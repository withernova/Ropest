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

        var entities = _softBodyQuery.ToEntityArray(Allocator.Temp);

        for (int e = 0; e < entities.Length; e++)
        {
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
            var resetDistHandle = resetDistJob.Schedule(edgeCount, 64, state.Dependency);

            var resetVolJob = new SoftBodyResetVolumeLambdaJob { Lambdas = volLambdaArr };
            var resetVolHandle = resetVolJob.Schedule(tetCount, 64, state.Dependency);

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

            state.Dependency = postSolveHandle;
        }

        entities.Dispose();
    }
}
