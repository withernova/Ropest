using Unity.Burst;
using Unity.Collections;
using Unity.Entities;
using Unity.Jobs;
using Unity.Mathematics;

/// <summary>
/// Cloth XPBD模拟系统 - 在FixedStep中运行
/// 调度顺序：ResetLambda -> PreSolve -> [SubSteps: Distance] -> PostSolve
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

        var entities = _clothQuery.ToEntityArray(Allocator.Temp);

        for (int e = 0; e < entities.Length; e++)
        {
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

            // 提取边数据到临时NativeArray
            int edgeCount = edges.Length;
            var edgeIndexA = new NativeArray<int>(edgeCount, Allocator.TempJob);
            var edgeIndexB = new NativeArray<int>(edgeCount, Allocator.TempJob);
            var restLengths = new NativeArray<float>(edgeCount, Allocator.TempJob);

            for (int i = 0; i < edgeCount; i++)
            {
                var edge = edges[i];
                edgeIndexA[i] = edge.IndexA;
                edgeIndexB[i] = edge.IndexB;
                restLengths[i] = edge.RestLength;
            }

            // === 1. 重置Lambda ===
            var resetJob = new ClothResetLambdaJob
            {
                Lambdas = lambdaArr
            };
            var resetHandle = resetJob.Schedule(edgeCount, 64, state.Dependency);

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

            // === 3. SubSteps: Distance Constraint + 自碰撞约束（空间哈希加速） ===
            // 创建空间哈希表（容量 = 粒子数 * 2 避免哈希冲突）
            var spatialHashMap = new NativeParallelMultiHashMap<int, int>(numParticles * 2, Allocator.TempJob);

            JobHandle constraintHandle = preSolveHandle;
            for (int step = 0; step < cfg.NumSubSteps; step++)
            {
                var distJob = new ClothDistanceConstraintJob
                {
                    Positions = posArr,
                    InvMasses = invMassArr,
                    EdgeIndexA = edgeIndexA,
                    EdgeIndexB = edgeIndexB,
                    RestLengths = restLengths,
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

                    // 第二步：并行自碰撞检测与修正
                    var selfCollisionJob = new ClothSelfCollisionJob
                    {
                        Positions = posArr,
                        InvMasses = invMassArr,
                        HashMap = spatialHashMap,
                        MinDistance = cfg.CollisionRadius,
                        CellSize = cellSize,
                        Subdivision = cfg.Subdivision + 1
                    };
                    constraintHandle = selfCollisionJob.Schedule(numParticles, 32, constraintHandle);
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

            // 释放临时数组
            edgeIndexA.Dispose(postSolveHandle);
            edgeIndexB.Dispose(postSolveHandle);
            restLengths.Dispose(postSolveHandle);
            spatialHashMap.Dispose(postSolveHandle);

            state.Dependency = postSolveHandle;
        }

        entities.Dispose();
    }
}
