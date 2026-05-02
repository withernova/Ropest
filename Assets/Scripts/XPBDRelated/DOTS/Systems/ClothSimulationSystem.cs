using Unity.Burst;
using Unity.Collections;
using Unity.Entities;
using Unity.Jobs;
using Unity.Mathematics;

[UpdateInGroup(typeof(FixedStepSimulationSystemGroup))]
public partial struct ClothSimulationSystem : ISystem
{
    private EntityQuery _clothQuery;

    public void OnCreate(ref SystemState state)
    {
        _clothQuery = state.GetEntityQuery(
            ComponentType.ReadOnly<ClothTag>(),
            ComponentType.ReadOnly<ClothSolverConfig>(),
            ComponentType.ReadWrite<ParticlePosition>(),
            ComponentType.ReadWrite<ParticlePrevPosition>(),
            ComponentType.ReadWrite<ParticleVelocity>(),
            ComponentType.ReadOnly<ParticleInvMass>(),
            ComponentType.ReadWrite<XPBDEdge>(),
            ComponentType.ReadWrite<XPBDDistanceLambda>()
        );
        state.RequireForUpdate(_clothQuery);
    }

    public void OnUpdate(ref SystemState state)
    {
        float dt = SystemAPI.Time.DeltaTime;
        if (dt <= 0f) return;

        _clothQuery.CompleteDependency();

        var entities = _clothQuery.ToEntityArray(Allocator.Temp);

        // 累积每个 cloth 实体 Job 链的依赖，最后统一写回 state.Dependency，
        // 避免循环中用 "state.Dependency = postSolveHandle" 覆盖掉前一个 cloth 的依赖。
        JobHandle combinedDependency = state.Dependency;

        for (int e = 0; e < entities.Length; e++)
        {
            combinedDependency.Complete();

            var entity = entities[e];
            var cfg = state.EntityManager.GetComponentData<ClothSolverConfig>(entity);
            var baseCfg = cfg.Base;
            int numParticles = baseCfg.NumParticles;

            var positions = state.EntityManager.GetBuffer<ParticlePosition>(entity);

            // 数据有效性检查：Buffer可能还未被填充
            if (positions.Length == 0 || numParticles <= 0) continue;

            var prevPositions = state.EntityManager.GetBuffer<ParticlePrevPosition>(entity);
            var velocities = state.EntityManager.GetBuffer<ParticleVelocity>(entity);
            var invMasses = state.EntityManager.GetBuffer<ParticleInvMass>(entity);
            var edges = state.EntityManager.GetBuffer<XPBDEdge>(entity);
            var lambdas = state.EntityManager.GetBuffer<XPBDDistanceLambda>(entity);

            var posArr = positions.Reinterpret<float3>().AsNativeArray();
            var prevArr = prevPositions.Reinterpret<float3>().AsNativeArray();
            var velArr = velocities.Reinterpret<float3>().AsNativeArray();
            var invMassArr = invMasses.Reinterpret<float>().AsNativeArray();
            var lambdaArr = lambdas.Reinterpret<float>().AsNativeArray();

            // 零拷贝获取 Edge Buffer（边数据是静态的，无需每帧重建NativeArray）
            var edgeArr = edges.AsNativeArray();
            int edgeCount = edges.Length;

            // === 1. 重置Lambda（共享 XPBDResetFloatBufferJob） ===
            var resetJob = new XPBDResetFloatBufferJob
            {
                Values = lambdaArr
            };
            var resetHandle = resetJob.Schedule(edgeCount, 64, combinedDependency);

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

            // === 3. SubSteps: Distance Constraint + 自碰撞 + 解析碰撞 ===
            // 创建空间哈希表（容量 = 粒子数 * 2 避免哈希冲突）
            var spatialHashMap = new NativeParallelMultiHashMap<int, int>(numParticles * 2, Allocator.TempJob);

            // 获取解析碰撞体数据（从全局管理器）
            var colliderData = AnalyticalColliderManager.GetColliderDataForJobs(Allocator.TempJob);

            JobHandle constraintHandle = preSolveHandle;
            for (int step = 0; step < baseCfg.NumSubSteps; step++)
            {
                // === 距离约束 ===
                // 根据开关在两种方案之间切换：
                //   - 图着色并行：按颜色组串行调度，每组内部用 IJobParallelFor 并行；
                //   - 串行 IJob：保持旧行为（全部边顺序求解）。
                if (cfg.UseGraphColoring && state.EntityManager.HasBuffer<XPBDEdgeColorRange>(entity))
                {
                    var colorRanges = state.EntityManager.GetBuffer<XPBDEdgeColorRange>(entity);
                    for (int c = 0; c < colorRanges.Length; c++)
                    {
                        int colorCount = colorRanges[c].Count;
                        if (colorCount <= 0) continue;

                        var distColoredJob = new XPBDDistanceConstraintColoredJob
                        {
                            Positions = posArr,
                            InvMasses = invMassArr,
                            Edges = edgeArr,
                            Lambdas = lambdaArr,
                            Stiffness = baseCfg.DistanceStiffness,
                            Dt = dt,
                            ColorStart = colorRanges[c].Start
                        };
                        constraintHandle = distColoredJob.Schedule(colorCount, 64, constraintHandle);
                    }
                }
                else
                {
                    var distJob = new XPBDDistanceConstraintJob
                    {
                        Positions = posArr,
                        InvMasses = invMassArr,
                        Edges = edgeArr,
                        Lambdas = lambdaArr,
                        Stiffness = baseCfg.DistanceStiffness,
                        Dt = dt
                    };
                    constraintHandle = distJob.Schedule(constraintHandle);
                }

                // 自碰撞：每隔2个SubStep做一次（平衡性能与效果） —— 布料专属
                if (step % 2 == 0)
                {
                    float cellSize = math.max(baseCfg.CollisionRadius * 2f, 0.01f);

                    // 第一步：构建空间哈希表（单线程）
                    var buildHashJob = new BuildSpatialHashJob
                    {
                        Positions = posArr,
                        HashMap = spatialHashMap,
                        CellSize = cellSize,
                        NumParticles = numParticles
                    };
                    constraintHandle = buildHashJob.Schedule(constraintHandle);

                    // 第二步：自碰撞检测与双向修正（IJob单线程，保证双向修正正确性）
                    var selfCollisionJob = new ClothSelfCollisionJob
                    {
                        Positions = posArr,
                        InvMasses = invMassArr,
                        HashMap = spatialHashMap,
                        MinDistance = baseCfg.CollisionRadius,
                        CellSize = cellSize,
                        Subdivision = cfg.Subdivision + 1,
                        NumParticles = numParticles
                    };
                    constraintHandle = selfCollisionJob.Schedule(constraintHandle);
                }

                // 解析碰撞：每个SubStep都做，与约束交替迭代（共享 XPBDAnalyticalCollisionJob）
                // 运行时开关：升起/落下阶段置 SkipAnalyticalCollision=true，避免布料穿过地面产生怪异网格
                if (colliderData.Length > 0 && !cfg.SkipAnalyticalCollision)
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
                        EnableFriction = false // 布料不启用 Job 内摩擦（保持旧行为一致）
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

            // 释放临时数组（edgeArr 是零拷贝视图无需释放）
            spatialHashMap.Dispose(postSolveHandle);
            colliderData.Dispose(postSolveHandle);

            // 不再直接覆盖 state.Dependency，而是累积到 combinedDependency。
            combinedDependency = JobHandle.CombineDependencies(combinedDependency, postSolveHandle);
        }

        state.Dependency = combinedDependency;

        entities.Dispose();
    }
}
