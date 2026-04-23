using Unity.Burst;
using Unity.Collections;
using Unity.Entities;
using Unity.Jobs;
using Unity.Mathematics;

/// <summary>
/// SoftBody XPBD模拟系统 - 在FixedStep中运行
/// 调度顺序：ResetLambda -> PreSolve -> [SubSteps: Distance + Volume + AnalyticalCollision] -> PostSolve
/// 可变形软体使用四面体体积约束 + 边距离约束来维持形状
/// 不需要自碰撞：距离约束+体积约束已足够防止自身穿模
/// </summary>
[UpdateInGroup(typeof(FixedStepSimulationSystemGroup))]
public partial struct SoftBodySimulationSystem : ISystem
{
    private EntityQuery _softBodyQuery;

    public void OnCreate(ref SystemState state)
    {
        _softBodyQuery = state.GetEntityQuery(
            ComponentType.ReadOnly<SoftBodyTag>(),
            ComponentType.ReadOnly<SoftBodySolverConfig>(),
            ComponentType.ReadWrite<ParticlePosition>(),
            ComponentType.ReadWrite<ParticlePrevPosition>(),
            ComponentType.ReadWrite<ParticleVelocity>(),
            ComponentType.ReadOnly<ParticleInvMass>(),
            ComponentType.ReadWrite<SoftBodyEdge>(),
            ComponentType.ReadWrite<SoftBodyDistanceLambda>(),
            ComponentType.ReadWrite<Tetrahedron>(),
            ComponentType.ReadWrite<TetrahedronRestVolume>(),
            ComponentType.ReadWrite<TetrahedronVolumeLambda>()
        );
        state.RequireForUpdate(_softBodyQuery);
    }

    public void OnUpdate(ref SystemState state)
    {
        float dt = SystemAPI.Time.DeltaTime;
        if (dt <= 0f) return;

        // 关键：在主线程上用 EntityManager.GetBuffer 访问 ParticlePosition 等组件之前，
        // 必须完成所有相关 ComponentType 上的未决读写 Job。
        // 不用 EntityManager.CompleteAllTrackedJobs()：它会 ClearDependencies()，
        // 破坏 ECS 的 fence 追踪，会波及其他系统（例如 AnalyticalCollider 相关 Job 链）。
        // 用 _softBodyQuery.CompleteDependency() 只等 Query 声明的 ComponentType 的 fence，精准无副作用。
        _softBodyQuery.CompleteDependency();

        var entities = _softBodyQuery.ToEntityArray(Allocator.Temp);

        // 累积每个 softbody 实体 Job 链的依赖，最后统一写回 state.Dependency，
        // 避免循环中直接赋值覆盖前一个实体的依赖。
        JobHandle combinedDependency = state.Dependency;

        for (int e = 0; e < entities.Length; e++)
        {
            // 关键：在用 EntityManager.GetBuffer 同步访问 ParticlePosition 之前，
            // 必须把上一轮循环里针对该 ComponentType 调度的 Job 全部完成。
            // ECS 的 AtomicSafety 按 ComponentType 做全局检查，不区分 entity。
            combinedDependency.Complete();

            var entity = entities[e];
            var cfg = state.EntityManager.GetComponentData<SoftBodySolverConfig>(entity);
            int numParticles = cfg.NumParticles;

            var positions = state.EntityManager.GetBuffer<ParticlePosition>(entity);

            // 数据有效性检查
            if (positions.Length == 0 || numParticles <= 0) continue;

            var prevPositions = state.EntityManager.GetBuffer<ParticlePrevPosition>(entity);
            var velocities = state.EntityManager.GetBuffer<ParticleVelocity>(entity);
            var invMasses = state.EntityManager.GetBuffer<ParticleInvMass>(entity);
            var edges = state.EntityManager.GetBuffer<SoftBodyEdge>(entity);
            var distLambdas = state.EntityManager.GetBuffer<SoftBodyDistanceLambda>(entity);
            var tets = state.EntityManager.GetBuffer<Tetrahedron>(entity);
            var tetRestVols = state.EntityManager.GetBuffer<TetrahedronRestVolume>(entity);
            var volLambdas = state.EntityManager.GetBuffer<TetrahedronVolumeLambda>(entity);

            var posArr = positions.Reinterpret<float3>().AsNativeArray();
            var prevArr = prevPositions.Reinterpret<float3>().AsNativeArray();
            var velArr = velocities.Reinterpret<float3>().AsNativeArray();
            var invMassArr = invMasses.Reinterpret<float>().AsNativeArray();
            var distLambdaArr = distLambdas.Reinterpret<float>().AsNativeArray();
            var volLambdaArr = volLambdas.Reinterpret<float>().AsNativeArray();

            // 零拷贝获取结构体Buffer（edges/tets/restVols是静态数据，无需每帧重建NativeArray）
            var edgeArr = edges.AsNativeArray();
            var tetArr = tets.AsNativeArray();
            var restVolArr = tetRestVols.AsNativeArray();

            int edgeCount = edges.Length;
            int tetCount = tets.Length;

            // === 1. 重置Lambda ===
            var resetDistJob = new SoftBodyResetDistanceLambdaJob { Lambdas = distLambdaArr };
            // 注意：依赖必须接在 combinedDependency 之后，避免多个 SoftBody 实体并行写相同 ComponentType。
            var resetDistHandle = resetDistJob.Schedule(edgeCount, 64, combinedDependency);

            var resetVolJob = new SoftBodyResetVolumeLambdaJob { Lambdas = volLambdaArr };
            var resetVolHandle = resetVolJob.Schedule(tetCount, 64, combinedDependency);

            var resetHandle = JobHandle.CombineDependencies(resetDistHandle, resetVolHandle);

            // === 2. PreSolve ===
            var preSolveJob = new SoftBodyPreSolveJob
            {
                Positions = posArr,
                PrevPositions = prevArr,
                Velocities = velArr,
                InvMasses = invMassArr,
                Gravity = cfg.Gravity,
                Dt = dt
            };
            var preSolveHandle = preSolveJob.Schedule(numParticles, 64, resetHandle);

            // === 3. SubSteps: Distance + Volume + 解析碰撞 ===
            var colliderData = AnalyticalColliderManager.GetColliderDataForJobs(Allocator.TempJob);

            JobHandle constraintHandle = preSolveHandle;
            for (int step = 0; step < cfg.NumSubSteps; step++)
            {
                // 距离约束
                var distJob = new SoftBodyDistanceConstraintJob
                {
                    Positions = posArr,
                    InvMasses = invMassArr,
                    Edges = edgeArr,
                    Lambdas = distLambdaArr,
                    Stiffness = cfg.DistanceStiffness,
                    Dt = dt
                };
                constraintHandle = distJob.Schedule(constraintHandle);

                // 体积约束
                var volJob = new SoftBodyVolumeConstraintJob
                {
                    Positions = posArr,
                    InvMasses = invMassArr,
                    Tets = tetArr,
                    RestVolumes = restVolArr,
                    Lambdas = volLambdaArr,
                    Stiffness = cfg.VolumeStiffness,
                    Dt = dt
                };
                constraintHandle = volJob.Schedule(constraintHandle);

                // 解析碰撞
                if (colliderData.Length > 0)
                {
                    var analyticalCollisionJob = new SoftBodyAnalyticalCollisionJob
                    {
                        Positions = posArr,
                        PrevPositions = prevArr,
                        InvMasses = invMassArr,
                        Colliders = colliderData,
                        ParticleRadius = cfg.CollisionRadius,
                        Friction = cfg.Friction,
                        NumParticles = numParticles
                    };
                    constraintHandle = analyticalCollisionJob.Schedule(constraintHandle);
                }
            }

            // === 4. PostSolve（含阻尼） ===
            var postSolveJob = new SoftBodyPostSolveJob
            {
                Positions = posArr,
                PrevPositions = prevArr,
                Velocities = velArr,
                InvMasses = invMassArr,
                OneOverDt = 1f / dt,
                Dt = dt,
                Damping = cfg.Damping
            };
            var postSolveHandle = postSolveJob.Schedule(numParticles, 64, constraintHandle);

            // 释放临时数组（仅 colliderData，其余都是零拷贝Buffer视图无需释放）
            colliderData.Dispose(postSolveHandle);

            combinedDependency = JobHandle.CombineDependencies(combinedDependency, postSolveHandle);
        }

        state.Dependency = combinedDependency;

        entities.Dispose();
    }
}
