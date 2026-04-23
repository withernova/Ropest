using Unity.Burst;
using Unity.Collections;
using Unity.Entities;
using Unity.Jobs;
using Unity.Mathematics;

/// <summary>
/// Rope XPBD模拟系统 - 在FixedStep中运行
/// 调度顺序：ResetLambda -> PreSolve -> [SubSteps: Edge + BendTwist] -> PostSolve
/// 注意：SystemAPI.Query最多支持8个类型参数，因此使用EntityQuery手动获取Buffer
/// </summary>
[UpdateInGroup(typeof(FixedStepSimulationSystemGroup))]
public partial struct RopeSimulationSystem : ISystem
{
    private EntityQuery _ropeQuery;

    public void OnCreate(ref SystemState state)
    {
        _ropeQuery = state.GetEntityQuery(
            ComponentType.ReadOnly<RopeTag>(),
            ComponentType.ReadOnly<RopeSolverConfig>(),
            ComponentType.ReadWrite<ParticlePosition>(),
            ComponentType.ReadWrite<ParticlePrevPosition>(),
            ComponentType.ReadWrite<ParticleVelocity>(),
            ComponentType.ReadOnly<ParticleInvMass>(),
            ComponentType.ReadWrite<GhostPosition>(),
            ComponentType.ReadWrite<GhostPrevPosition>(),
            ComponentType.ReadWrite<GhostVelocity>(),
            ComponentType.ReadOnly<GhostInvMass>(),
            ComponentType.ReadOnly<RopeRestLength>(),
            ComponentType.ReadOnly<InitDarbouxVector>(),
            ComponentType.ReadWrite<EdgeLambda0>(),
            ComponentType.ReadWrite<EdgeLambda1>(),
            ComponentType.ReadWrite<EdgeLambda2>(),
            ComponentType.ReadWrite<BendTwistLambda>()
        );
        state.RequireForUpdate(_ropeQuery);
    }

    public void OnUpdate(ref SystemState state)
    {
        float dt = SystemAPI.Time.DeltaTime;
        if (dt <= 0f) return;

        // 拉取场景中的解析碰撞体（与 Cloth/SoftBody 相同的数据源）。
        // 每帧拷一份供所有 Rope 共用，PostSolve 后再 Dispose。
        var colliderData = AnalyticalColliderManager.GetColliderDataForJobs(Allocator.TempJob);

        var entities = _ropeQuery.ToEntityArray(Allocator.Temp);

        for (int e = 0; e < entities.Length; e++)
        {
            var entity = entities[e];
            var cfg = state.EntityManager.GetComponentData<RopeSolverConfig>(entity);
            int numPoints = cfg.NumPoints;
            int numGhosts = cfg.NumGhostPoints;

            // 获取所有DynamicBuffer
            var positions = state.EntityManager.GetBuffer<ParticlePosition>(entity);

            // 数据有效性检查：Buffer可能还未被填充
            if (positions.Length == 0 || numPoints <= 0) continue;
            var prevPositions = state.EntityManager.GetBuffer<ParticlePrevPosition>(entity);
            var velocities = state.EntityManager.GetBuffer<ParticleVelocity>(entity);
            var invMasses = state.EntityManager.GetBuffer<ParticleInvMass>(entity);
            var ghostPositions = state.EntityManager.GetBuffer<GhostPosition>(entity);
            var ghostPrevPositions = state.EntityManager.GetBuffer<GhostPrevPosition>(entity);
            var ghostVelocities = state.EntityManager.GetBuffer<GhostVelocity>(entity);
            var ghostInvMasses = state.EntityManager.GetBuffer<GhostInvMass>(entity);
            var restLengths = state.EntityManager.GetBuffer<RopeRestLength>(entity);
            var initDarboux = state.EntityManager.GetBuffer<InitDarbouxVector>(entity);
            var edgeLambda0 = state.EntityManager.GetBuffer<EdgeLambda0>(entity);
            var edgeLambda1 = state.EntityManager.GetBuffer<EdgeLambda1>(entity);
            var edgeLambda2 = state.EntityManager.GetBuffer<EdgeLambda2>(entity);
            var bendLambdas = state.EntityManager.GetBuffer<BendTwistLambda>(entity);

            // 获取NativeArray引用
            var posArr = positions.Reinterpret<float3>().AsNativeArray();
            var prevArr = prevPositions.Reinterpret<float3>().AsNativeArray();
            var velArr = velocities.Reinterpret<float3>().AsNativeArray();
            var invMassArr = invMasses.Reinterpret<float>().AsNativeArray();
            var ghostPosArr = ghostPositions.Reinterpret<float3>().AsNativeArray();
            var ghostPrevArr = ghostPrevPositions.Reinterpret<float3>().AsNativeArray();
            var ghostVelArr = ghostVelocities.Reinterpret<float3>().AsNativeArray();
            var ghostInvMassArr = ghostInvMasses.Reinterpret<float>().AsNativeArray();
            var restLenArr = restLengths.Reinterpret<float>().AsNativeArray();
            var initDarbouxArr = initDarboux.Reinterpret<float3>().AsNativeArray();
            var el0 = edgeLambda0.Reinterpret<float>().AsNativeArray();
            var el1 = edgeLambda1.Reinterpret<float>().AsNativeArray();
            var el2 = edgeLambda2.Reinterpret<float>().AsNativeArray();
            var btLambdas = bendLambdas.Reinterpret<float3>().AsNativeArray();

            // === 1. 重置Lambda ===
            var resetEdgeJob = new RopeResetLambdaJob
            {
                EdgeLambda0 = el0,
                EdgeLambda1 = el1,
                EdgeLambda2 = el2
            };
            var resetEdgeHandle = resetEdgeJob.Schedule(numPoints, 64, state.Dependency);

            int numBend = math.max(0, numPoints - 2);
            var resetBendJob = new RopeResetBendLambdaJob
            {
                BendLambdas = btLambdas
            };
            var resetBendHandle = resetBendJob.Schedule(numBend, 64, state.Dependency);

            var resetHandle = JobHandle.CombineDependencies(resetEdgeHandle, resetBendHandle);

            // === 2. PreSolve（拆分为Point和Ghost两个并行Job） ===
            var preSolvePointJob = new RopePreSolvePointJob
            {
                PointPos = posArr,
                PrevPos = prevArr,
                Vel = velArr,
                PointInvMass = invMassArr,
                Gravity = cfg.Gravity,
                GravityFactor = cfg.GravityFactor,
                Dt = dt
            };
            var preSolvePointHandle = preSolvePointJob.Schedule(numPoints, 64, resetHandle);

            var preSolveGhostJob = new RopePreSolveGhostJob
            {
                GhostPos = ghostPosArr,
                GhostPrev = ghostPrevArr,
                GhostVels = ghostVelArr,
                Dt = dt
            };
            var preSolveGhostHandle = preSolveGhostJob.Schedule(numGhosts, 64, resetHandle);

            var preSolveHandle = JobHandle.CombineDependencies(preSolvePointHandle, preSolveGhostHandle);

            // === 3. SubSteps: Edge + BendTwist ===
            JobHandle constraintHandle = preSolveHandle;
            for (int step = 0; step < cfg.NumSubSteps; step++)
            {
                var edgeJob = new RopeEdgeConstraintJob
                {
                    PointPos = posArr,
                    GhostPos = ghostPosArr,
                    PointInvMass = invMassArr,
                    GhostInvMass = ghostInvMassArr,
                    RestLengths = restLenArr,
                    Lambda0 = el0,
                    Lambda1 = el1,
                    Lambda2 = el2,
                    Stiffness = cfg.EdgeStiffness,
                    GhostDistance = cfg.GhostDistance,
                    Dt = dt,
                    NumPoints = numPoints
                };
                constraintHandle = edgeJob.Schedule(constraintHandle);

                if (numBend > 0)
                {
                    var bendJob = new RopeBendTwistConstraintJob
                    {
                        PointPos = posArr,
                        GhostPos = ghostPosArr,
                        PointInvMass = invMassArr,
                        GhostInvMass = ghostInvMassArr,
                        InitDarboux = initDarbouxArr,
                        RestLengths = restLenArr,
                        Lambdas = btLambdas,
                        BendTwistKs = cfg.BendTwistKs,
                        // 修复：BendTwist 必须有自己独立的 compliance（原版 RopeContraints.cs 中基类默认 stiff=0.6f）。
                        // 先前错误地共用了 cfg.EdgeStiffness=0.0001f，导致 α 比原版小 ~6000 倍，
                        // 约束过于刚性 → factor_matrix 数值病态 → 第2~4帧发散为 NaN/Inf，绳子消失。
                        Stiffness = cfg.BendTwistStiffness,
                        Dt = dt,
                        NumPoints = numPoints
                    };
                    constraintHandle = bendJob.Schedule(constraintHandle);
                }

                // === 3c. 解析碰撞（每个 substep 内解一次，避免碰撞修正破坏 Edge 约束） ===
                if (colliderData.Length > 0)
                {
                    var collisionJob = new RopeAnalyticalCollisionJob
                    {
                        PointPos = posArr,
                        PrevPos = prevArr,
                        PointInvMass = invMassArr,
                        Colliders = colliderData,
                        ParticleRadius = cfg.Radius,
                        Friction = cfg.Friction,
                        NumPoints = numPoints
                    };
                    constraintHandle = collisionJob.Schedule(constraintHandle);
                }
            }

            // === 4. PostSolve（拆分为两步：先更新Vel+阻尼，再更新GhostVels） ===
            var postSolveVelJob = new RopePostSolveVelocityJob
            {
                PointPos = posArr,
                PrevPos = prevArr,
                Vel = velArr,
                OneOverDt = 1f / dt,
                Dt = dt,
                Damping = cfg.Damping
            };
            var velHandle = postSolveVelJob.Schedule(numPoints, 64, constraintHandle);

            var postSolveGhostVelJob = new RopePostSolveGhostVelJob
            {
                Vel = velArr,
                GhostVels = ghostVelArr
            };
            state.Dependency = postSolveGhostVelJob.Schedule(numGhosts, 64, velHandle);
        }

        entities.Dispose();

        // colliderData 在所有 Rope 的 Job 完成后再释放；
        // 由于每次 entity 循环都把依赖链挂到 state.Dependency，这里 Dispose(state.Dependency) 即可安全释放。
        colliderData.Dispose(state.Dependency);
    }
}
