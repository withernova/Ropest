using Unity.Collections;
using Unity.Entities;
using Unity.Jobs;
using Unity.Mathematics;
using UnityEngine;

/// <summary>
/// 跨体碰撞系统 - 处理不同实体（Cloth/SoftBody）之间的粒子-粒子碰撞
/// 调度位置：在 ClothSimulationSystem 和 SoftBodySimulationSystem 之后、XPBDCollisionSystem 之前
///
/// 工作流程（每帧）：
///   1) 主线程：收集所有 Cloth + SoftBody 实体的 Buffer 数据到全局扁平 NativeArray
///   2) 迭代 N 次：BuildHash → Resolve → Apply
///   3) 主线程：把"净修正 = GlobalPositions - OriginalPositions"写回各实体 Buffer
///      同时同步更新 Velocity = (new_pos - prev_pos) / dt，让修正立刻体现在速度上
///      （否则下一帧 PreSolve 基于旧速度积分，视觉上会"延迟"甚至"不动"）
///
/// 设计要点：
///   - 所有 Buffer 访问都在主线程完成，避免"边 Schedule 边 GetBuffer"的 AtomicSafety 冲突
///   - Velocity 同步更新是让双向形变耦合立刻可见的关键
/// </summary>
[UpdateInGroup(typeof(FixedStepSimulationSystemGroup))]
[UpdateAfter(typeof(ClothSimulationSystem))]
[UpdateAfter(typeof(SoftBodySimulationSystem))]
[UpdateBefore(typeof(XPBDCollisionSystem))]
public partial struct CrossBodyCollisionSystem : ISystem
{
    private EntityQuery _clothQuery;
    private EntityQuery _softBodyQuery;

    // 迭代次数：越高越稳定但越慢；2 次是性能和稳定性的良好平衡
    private const int NumIterations = 2;
    // 跨体摩擦：用于耗散切向速度，防止粒子在对方表面无限滑动
    private const float CrossBodyFriction = 0.2f;

    // 诊断：每 N 帧打一次日志，观察跨体碰撞是否真的在工作、修正幅度多大
    // 改成 0 可关闭日志
    private const int DebugLogInterval = 60;

    public void OnCreate(ref SystemState state)
    {
        _clothQuery = state.GetEntityQuery(
            ComponentType.ReadOnly<ClothTag>(),
            ComponentType.ReadOnly<ClothSolverConfig>(),
            ComponentType.ReadWrite<ParticlePosition>(),
            ComponentType.ReadWrite<ParticlePrevPosition>(),
            ComponentType.ReadWrite<ParticleVelocity>(),
            ComponentType.ReadOnly<ParticleInvMass>()
        );
        _softBodyQuery = state.GetEntityQuery(
            ComponentType.ReadOnly<SoftBodyTag>(),
            ComponentType.ReadOnly<SoftBodySolverConfig>(),
            ComponentType.ReadWrite<ParticlePosition>(),
            ComponentType.ReadWrite<ParticlePrevPosition>(),
            ComponentType.ReadWrite<ParticleVelocity>(),
            ComponentType.ReadOnly<ParticleInvMass>()
        );
    }

    public void OnUpdate(ref SystemState state)
    {
        float dt = SystemAPI.Time.DeltaTime;
        if (dt <= 0f) return;

        var clothEntities = _clothQuery.ToEntityArray(Allocator.Temp);
        var softBodyEntities = _softBodyQuery.ToEntityArray(Allocator.Temp);

        int totalBodies = clothEntities.Length + softBodyEntities.Length;

        if (totalBodies < 2)
        {
            clothEntities.Dispose();
            softBodyEntities.Dispose();
            return;
        }

        // 关键：在主线程访问 Buffer 前，先完成所有相关 ComponentType 上的 Job
        state.CompleteDependency();
        _clothQuery.CompleteDependency();
        _softBodyQuery.CompleteDependency();

        int numCloth = clothEntities.Length;
        int numSoft = softBodyEntities.Length;

        var offsets = new NativeArray<int>(totalBodies, Allocator.Temp);
        var sizes = new NativeArray<int>(totalBodies, Allocator.Temp);
        var radii = new NativeArray<float>(totalBodies, Allocator.Temp);
        var frictions = new NativeArray<float>(totalBodies, Allocator.Temp);

        int cursor = 0;
        float maxRadius = 0f;

        for (int i = 0; i < numCloth; i++)
        {
            var cfg = state.EntityManager.GetComponentData<ClothSolverConfig>(clothEntities[i]);
            offsets[i] = cursor;
            sizes[i] = cfg.NumParticles;
            radii[i] = cfg.CollisionRadius;
            frictions[i] = math.max(cfg.Friction, CrossBodyFriction);
            cursor += cfg.NumParticles;
            if (cfg.CollisionRadius > maxRadius) maxRadius = cfg.CollisionRadius;
        }
        for (int i = 0; i < numSoft; i++)
        {
            var cfg = state.EntityManager.GetComponentData<SoftBodySolverConfig>(softBodyEntities[i]);
            int idx = numCloth + i;
            offsets[idx] = cursor;
            sizes[idx] = cfg.NumParticles;
            radii[idx] = cfg.CollisionRadius;
            frictions[idx] = math.max(cfg.Friction, CrossBodyFriction);
            cursor += cfg.NumParticles;
            if (cfg.CollisionRadius > maxRadius) maxRadius = cfg.CollisionRadius;
        }

        int totalParticles = cursor;

        if (totalParticles <= 0 || maxRadius <= 0f)
        {
            clothEntities.Dispose();
            softBodyEntities.Dispose();
            offsets.Dispose();
            sizes.Dispose();
            radii.Dispose();
            frictions.Dispose();
            return;
        }

        // === 分配全局扁平数组 ===
        var globalPositions = new NativeArray<float3>(totalParticles, Allocator.TempJob, NativeArrayOptions.UninitializedMemory);
        var globalInvMasses = new NativeArray<float>(totalParticles, Allocator.TempJob, NativeArrayOptions.UninitializedMemory);
        var globalRadii = new NativeArray<float>(totalParticles, Allocator.TempJob, NativeArrayOptions.UninitializedMemory);
        var globalBodyIds = new NativeArray<int>(totalParticles, Allocator.TempJob, NativeArrayOptions.UninitializedMemory);
        var globalCorrections = new NativeArray<float3>(totalParticles, Allocator.TempJob, NativeArrayOptions.ClearMemory);

        // === 1. 主线程收集：把所有实体的 Buffer 数据拷到全局数组 ===
        for (int i = 0; i < numCloth; i++)
        {
            var entity = clothEntities[i];
            var posBuf = state.EntityManager.GetBuffer<ParticlePosition>(entity, true);
            var invMassBuf = state.EntityManager.GetBuffer<ParticleInvMass>(entity, true);
            int num = sizes[i];
            int off = offsets[i];
            float r = radii[i];
            int bid = i;
            for (int k = 0; k < num; k++)
            {
                globalPositions[off + k] = posBuf[k].Value;
                globalInvMasses[off + k] = invMassBuf[k].Value;
                globalRadii[off + k] = r;
                globalBodyIds[off + k] = bid;
            }
        }
        for (int i = 0; i < numSoft; i++)
        {
            int idx = numCloth + i;
            var entity = softBodyEntities[i];
            var posBuf = state.EntityManager.GetBuffer<ParticlePosition>(entity, true);
            var invMassBuf = state.EntityManager.GetBuffer<ParticleInvMass>(entity, true);
            int num = sizes[idx];
            int off = offsets[idx];
            float r = radii[idx];
            int bid = idx;
            for (int k = 0; k < num; k++)
            {
                globalPositions[off + k] = posBuf[k].Value;
                globalInvMasses[off + k] = invMassBuf[k].Value;
                globalRadii[off + k] = r;
                globalBodyIds[off + k] = bid;
            }
        }

        // 原始位置快照
        var originalPositions = new NativeArray<float3>(globalPositions, Allocator.TempJob);

        float cellSize = math.max(maxRadius * 2f, 0.01f);
        var hashMap = new NativeParallelMultiHashMap<int, int>(totalParticles * 2, Allocator.TempJob);

        // === 2. 迭代：BuildHash → Resolve → Apply ===
        JobHandle iterHandle = default;
        for (int iter = 0; iter < NumIterations; iter++)
        {
            var buildHashJob = new BuildCrossBodyHashJob
            {
                GlobalPositions = globalPositions,
                HashMap = hashMap,
                CellSize = cellSize,
                NumParticles = totalParticles
            };
            iterHandle = buildHashJob.Schedule(iterHandle);

            var resolveJob = new CrossBodyCollisionResolveJob
            {
                GlobalPositions = globalPositions,
                GlobalInvMasses = globalInvMasses,
                GlobalRadii = globalRadii,
                GlobalBodyIds = globalBodyIds,
                HashMap = hashMap,
                GlobalCorrections = globalCorrections,
                CellSize = cellSize,
                NumParticles = totalParticles
            };
            iterHandle = resolveJob.Schedule(iterHandle);

            var applyJob = new ApplyGlobalCorrectionJob
            {
                GlobalPositions = globalPositions,
                GlobalCorrections = globalCorrections,
                NumParticles = totalParticles
            };
            iterHandle = applyJob.Schedule(iterHandle);
        }

        iterHandle.Complete();

        // === 3. 诊断（可选）：统计本帧跨体碰撞的净修正 ===
        int clothContacts = 0, softContacts = 0;
        float clothMaxCorr = 0f, softMaxCorr = 0f;
        bool doDebug = DebugLogInterval > 0 && (Time.frameCount % DebugLogInterval == 0);

        // === 4. 主线程写回：净修正 + 速度同步 + 摩擦 ===
        for (int i = 0; i < numCloth; i++)
        {
            var entity = clothEntities[i];
            var posBuf = state.EntityManager.GetBuffer<ParticlePosition>(entity);
            var prevBuf = state.EntityManager.GetBuffer<ParticlePrevPosition>(entity);
            var velBuf = state.EntityManager.GetBuffer<ParticleVelocity>(entity);
            int num = sizes[i];
            int off = offsets[i];
            float friction = frictions[i];
            ScatterBackToBuffer(posBuf, prevBuf, velBuf, globalPositions, originalPositions, globalInvMasses,
                off, num, friction, dt, doDebug, ref clothContacts, ref clothMaxCorr);
        }
        for (int i = 0; i < numSoft; i++)
        {
            int idx = numCloth + i;
            var entity = softBodyEntities[i];
            var posBuf = state.EntityManager.GetBuffer<ParticlePosition>(entity);
            var prevBuf = state.EntityManager.GetBuffer<ParticlePrevPosition>(entity);
            var velBuf = state.EntityManager.GetBuffer<ParticleVelocity>(entity);
            int num = sizes[idx];
            int off = offsets[idx];
            float friction = frictions[idx];
            ScatterBackToBuffer(posBuf, prevBuf, velBuf, globalPositions, originalPositions, globalInvMasses,
                off, num, friction, dt, doDebug, ref softContacts, ref softMaxCorr);
        }

        if (doDebug)
        {
            Debug.Log($"[CrossBodyCollision] totalParticles={totalParticles} cellSize={cellSize:F3} maxRadius={maxRadius:F3} | " +
                      $"cloth接触粒子={clothContacts} 最大修正={clothMaxCorr:F4} | " +
                      $"soft接触粒子={softContacts} 最大修正={softMaxCorr:F4}");
        }

        // === 清理 ===
        globalPositions.Dispose();
        originalPositions.Dispose();
        globalInvMasses.Dispose();
        globalRadii.Dispose();
        globalBodyIds.Dispose();
        globalCorrections.Dispose();
        hashMap.Dispose();

        clothEntities.Dispose();
        softBodyEntities.Dispose();
        offsets.Dispose();
        sizes.Dispose();
        radii.Dispose();
        frictions.Dispose();
    }

    /// <summary>
    /// 主线程：把净修正写回实体 Buffer，并同步更新速度（关键：让碰撞效果立即可见，而非延迟一帧）
    /// </summary>
    private static void ScatterBackToBuffer(
        DynamicBuffer<ParticlePosition> posBuf,
        DynamicBuffer<ParticlePrevPosition> prevBuf,
        DynamicBuffer<ParticleVelocity> velBuf,
        NativeArray<float3> globalPositions,
        NativeArray<float3> originalPositions,
        NativeArray<float> globalInvMasses,
        int offset, int num, float friction, float dt,
        bool countStats, ref int contactCount, ref float maxCorrLen)
    {
        float invDt = 1f / dt;
        for (int k = 0; k < num; k++)
        {
            int gi = offset + k;
            if (globalInvMasses[gi] == 0f) continue;

            float3 corr = globalPositions[gi] - originalPositions[gi];
            float corrLen = math.length(corr);
            if (corrLen < 1e-8f) continue;

            if (countStats)
            {
                contactCount++;
                if (corrLen > maxCorrLen) maxCorrLen = corrLen;
            }

            // 应用位置修正
            float3 newPos = posBuf[k].Value + corr;
            posBuf[k] = new ParticlePosition { Value = newPos };

            // 摩擦：沿切向把 prev 朝 pos 拉近
            float3 prev = prevBuf[k].Value;
            if (friction > 0f)
            {
                float3 normal = corr / corrLen;
                float3 disp = newPos - prev;
                float3 tangent = disp - math.dot(disp, normal) * normal;
                float tangentLen = math.length(tangent);

                if (tangentLen < 1e-4f)
                {
                    prev = newPos;
                }
                else
                {
                    float fric = math.saturate(friction);
                    prev += tangent * fric;
                }
                prevBuf[k] = new ParticlePrevPosition { Value = prev };
            }

            // 关键：同步更新速度为"新位置 - prev"/ dt
            // 这样下一帧 PreSolve 会基于修正后的速度积分，修正效果立即可见
            velBuf[k] = new ParticleVelocity { Value = (newPos - prev) * invDt };
        }
    }
}

/// <summary>
/// 把本次迭代累积的 Corrections 应用到 GlobalPositions，然后清零 Corrections
/// </summary>
[Unity.Burst.BurstCompile]
public struct ApplyGlobalCorrectionJob : Unity.Jobs.IJob
{
    public NativeArray<float3> GlobalPositions;
    public NativeArray<float3> GlobalCorrections;
    public int NumParticles;

    public void Execute()
    {
        for (int i = 0; i < NumParticles; i++)
        {
            GlobalPositions[i] += GlobalCorrections[i];
            GlobalCorrections[i] = float3.zero;
        }
    }
}