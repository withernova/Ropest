using Unity.Burst;
using Unity.Collections;
using Unity.Entities;
using Unity.Jobs;
using Unity.Mathematics;

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
            ComponentType.ReadWrite<XPBDEdge>(),
            ComponentType.ReadWrite<XPBDDistanceLambda>(),
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

        _softBodyQuery.CompleteDependency();

        var entities = _softBodyQuery.ToEntityArray(Allocator.Temp);

        // 累积每个 softbody 实体 Job 链的依赖，最后统一写回 state.Dependency，
        // 避免循环中直接赋值覆盖前一个实体的依赖。
        JobHandle combinedDependency = state.Dependency;

        for (int e = 0; e < entities.Length; e++)
        {
            combinedDependency.Complete();

            var entity = entities[e];
            var cfg = state.EntityManager.GetComponentData<SoftBodySolverConfig>(entity);
            var baseCfg = cfg.Base;
            int numParticles = baseCfg.NumParticles;

            var positions = state.EntityManager.GetBuffer<ParticlePosition>(entity);

            // 数据有效性检查
            if (positions.Length == 0 || numParticles <= 0) continue;

            var prevPositions = state.EntityManager.GetBuffer<ParticlePrevPosition>(entity);
            var velocities = state.EntityManager.GetBuffer<ParticleVelocity>(entity);
            var invMasses = state.EntityManager.GetBuffer<ParticleInvMass>(entity);
            var edges = state.EntityManager.GetBuffer<XPBDEdge>(entity);
            var distLambdas = state.EntityManager.GetBuffer<XPBDDistanceLambda>(entity);
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

            // === 1. 重置Lambda（共享 XPBDResetFloatBufferJob） ===
            var resetDistJob = new XPBDResetFloatBufferJob { Values = distLambdaArr };
            // 注意：依赖必须接在 combinedDependency 之后，避免多个 SoftBody 实体并行写相同 ComponentType。
            var resetDistHandle = resetDistJob.Schedule(edgeCount, 64, combinedDependency);

            var resetVolJob = new XPBDResetFloatBufferJob { Values = volLambdaArr };
            var resetVolHandle = resetVolJob.Schedule(tetCount, 64, combinedDependency);

            var resetHandle = JobHandle.CombineDependencies(resetDistHandle, resetVolHandle);

            // === 2. PreSolve（共享 XPBDPreSolveJob） ===
            var preSolveJob = new XPBDPreSolveJob
            {
                Positions = posArr,
                PrevPositions = prevArr,
                Velocities = velArr,
                InvMasses = invMassArr,
                Gravity = baseCfg.Gravity,
                Dt = dt
            };
            var preSolveHandle = preSolveJob.Schedule(numParticles, 64, resetHandle);

            // === 3. SubSteps: Distance + Volume + 解析碰撞 ===
            var colliderData = AnalyticalColliderManager.GetColliderDataForJobs(Allocator.TempJob);

            JobHandle constraintHandle = preSolveHandle;
            for (int step = 0; step < baseCfg.NumSubSteps; step++)
            {
                // 距离约束（共享 XPBDDistanceConstraintJob）
                var distJob = new XPBDDistanceConstraintJob
                {
                    Positions = posArr,
                    InvMasses = invMassArr,
                    Edges = edgeArr,
                    Lambdas = distLambdaArr,
                    Stiffness = baseCfg.DistanceStiffness,
                    Dt = dt
                };
                constraintHandle = distJob.Schedule(constraintHandle);

                // 体积约束（软体专属）
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

                // 解析碰撞（共享 XPBDAnalyticalCollisionJob，软体启用摩擦）
                if (colliderData.Length > 0)
                {
                    var analyticalCollisionJob = new XPBDAnalyticalCollisionJob
                    {
                        Positions = posArr,
                        PrevPositions = prevArr,
                        InvMasses = invMassArr,
                        Colliders = colliderData,
                        ParticleRadius = baseCfg.CollisionRadius,
                        Friction = baseCfg.Friction,
                        NumParticles = numParticles,
                        EnableFriction = true
                    };
                    constraintHandle = analyticalCollisionJob.Schedule(constraintHandle);
                }
            }

            // === 4. PostSolve（含阻尼，共享 XPBDPostSolveJob） ===
            var postSolveJob = new XPBDPostSolveJob
            {
                Positions = posArr,
                PrevPositions = prevArr,
                Velocities = velArr,
                InvMasses = invMassArr,
                OneOverDt = 1f / dt,
                Dt = dt,
                Damping = baseCfg.Damping
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
