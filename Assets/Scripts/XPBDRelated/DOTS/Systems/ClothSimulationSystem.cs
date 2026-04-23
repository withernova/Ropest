using Unity.Burst;
using Unity.Collections;
using Unity.Entities;
using Unity.Jobs;
using Unity.Mathematics;

/// <summary>
/// Cloth XPBD模拟系统 - 在FixedStep中运行
/// 调度顺序：ResetLambda -> PreSolve -> [SubSteps: Distance + SelfCollision + AnalyticalCollision] -> PostSolve
/// 场景碰撞使用解析碰撞（球体/Box），在SubStep内与约束交替迭代
/// </summary>
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
            ComponentType.ReadWrite<ClothEdge>(),
            ComponentType.ReadWrite<ClothDistanceLambda>()
        );
        state.RequireForUpdate(_clothQuery);
    }

    public void OnUpdate(ref SystemState state)
    {
        float dt = SystemAPI.Time.DeltaTime;
        if (dt <= 0f) return;

        // 关键：在主线程上用 EntityManager.GetBuffer 访问 ParticlePosition 等组件之前，
        // 必须完成所有相关 ComponentType 上的未决读写 Job。
        // 这里不能用 state.Dependency.Complete()——它只等本系统 Dependency 链里记录的依赖，
        // 无法等到"循环内通过 NativeArray 间接写 Buffer 的 Job"的 ComponentType 级 fence。
        // 也不用 EntityManager.CompleteAllTrackedJobs()——它会 ClearDependencies()，
        // 破坏 ECS 的 fence 追踪，导致后续系统无法正确同步（从而影响碰撞等交互）。
        // 正确做法：用 EntityQuery.CompleteDependency()，它只等 Query 声明的 ComponentType 上的 fence，
        // 精准且无副作用。
        _clothQuery.CompleteDependency();

        var entities = _clothQuery.ToEntityArray(Allocator.Temp);

        // 累积每个 cloth 实体 Job 链的依赖，最后统一写回 state.Dependency，
        // 避免循环中用 "state.Dependency = postSolveHandle" 覆盖掉前一个 cloth 的依赖。
        JobHandle combinedDependency = state.Dependency;

        for (int e = 0; e < entities.Length; e++)
        {
            // 关键：在用 EntityManager.GetBuffer 同步访问 ParticlePosition 之前，
            // 必须把上一轮循环里针对该 ComponentType 调度的 Job 全部完成。
            // ECS 的 AtomicSafety 按 ComponentType 做全局检查，不区分 entity，
            // 因此即使这里读的是不同 entity 的 buffer，只要有未完成的同类型写 Job 就会抛异常。
            combinedDependency.Complete();

            var entity = entities[e];
            var cfg = state.EntityManager.GetComponentData<ClothSolverConfig>(entity);
            int numParticles = cfg.NumParticles;

            var positions = state.EntityManager.GetBuffer<ParticlePosition>(entity);

            // 数据有效性检查：Buffer可能还未被填充
            if (positions.Length == 0 || numParticles <= 0) continue;

            var prevPositions = state.EntityManager.GetBuffer<ParticlePrevPosition>(entity);
            var velocities = state.EntityManager.GetBuffer<ParticleVelocity>(entity);
            var invMasses = state.EntityManager.GetBuffer<ParticleInvMass>(entity);
            var edges = state.EntityManager.GetBuffer<ClothEdge>(entity);
            var lambdas = state.EntityManager.GetBuffer<ClothDistanceLambda>(entity);

            var posArr = positions.Reinterpret<float3>().AsNativeArray();
            var prevArr = prevPositions.Reinterpret<float3>().AsNativeArray();
            var velArr = velocities.Reinterpret<float3>().AsNativeArray();
            var invMassArr = invMasses.Reinterpret<float>().AsNativeArray();
            var lambdaArr = lambdas.Reinterpret<float>().AsNativeArray();

            // 零拷贝获取 ClothEdge Buffer（边数据是静态的，无需每帧重建NativeArray）
            var edgeArr = edges.AsNativeArray();
            int edgeCount = edges.Length;

            // === 1. 重置Lambda ===
            var resetJob = new ClothResetLambdaJob
            {
                Lambdas = lambdaArr
            };
            // 注意：依赖必须接在 combinedDependency 之后，而不是 state.Dependency。
            // 因为 ParticlePosition 等 ComponentType 是全局共享的，上一轮循环里
            // 另一个 Cloth 实体调度的 Job 也在写此类型，ECS 不允许并行写，
            // 必须让多个 Cloth 形成串行链。
            var resetHandle = resetJob.Schedule(edgeCount, 64, combinedDependency);

            // === 2. PreSolve ===
            var preSolveJob = new ClothPreSolveJob
            {
                Positions = posArr,
                PrevPositions = prevArr,
                Velocities = velArr,
                InvMasses = invMassArr,
                Gravity = cfg.Gravity,
                Dt = dt
            };
            var preSolveHandle = preSolveJob.Schedule(numParticles, 64, resetHandle);

            // === 3. SubSteps: Distance Constraint + 自碰撞 + 解析碰撞 ===
            // 创建空间哈希表（容量 = 粒子数 * 2 避免哈希冲突）
            var spatialHashMap = new NativeParallelMultiHashMap<int, int>(numParticles * 2, Allocator.TempJob);

            // 获取解析碰撞体数据（从全局管理器）
            var colliderData = AnalyticalColliderManager.GetColliderDataForJobs(Allocator.TempJob);

            JobHandle constraintHandle = preSolveHandle;
            for (int step = 0; step < cfg.NumSubSteps; step++)
            {
                var distJob = new ClothDistanceConstraintJob
                {
                    Positions = posArr,
                    InvMasses = invMassArr,
                    Edges = edgeArr,
                    Lambdas = lambdaArr,
                    Stiffness = cfg.DistanceStiffness,
                    Dt = dt
                };
                constraintHandle = distJob.Schedule(constraintHandle);

                // 自碰撞：每隔2个SubStep做一次（平衡性能与效果）
                if (step % 2 == 0)
                {
                    float cellSize = math.max(cfg.CollisionRadius * 2f, 0.01f);

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
                        MinDistance = cfg.CollisionRadius,
                        CellSize = cellSize,
                        Subdivision = cfg.Subdivision + 1,
                        NumParticles = numParticles
                    };
                    constraintHandle = selfCollisionJob.Schedule(constraintHandle);
                }

                // 解析碰撞：每个SubStep都做，与约束交替迭代
                if (colliderData.Length > 0)
                {
                    var analyticalCollisionJob = new ClothAnalyticalCollisionJob
                    {
                        Positions = posArr,
                        InvMasses = invMassArr,
                        Colliders = colliderData,
                        ParticleRadius = cfg.CollisionRadius,
                        NumParticles = numParticles
                    };
                    constraintHandle = analyticalCollisionJob.Schedule(constraintHandle);
                }
            }

            // === 4. PostSolve（含阻尼） ===
            var postSolveJob = new ClothPostSolveJob
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
