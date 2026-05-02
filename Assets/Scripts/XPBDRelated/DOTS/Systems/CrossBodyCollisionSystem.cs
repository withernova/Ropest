using Unity.Collections;
using Unity.Entities;
using Unity.Jobs;
using Unity.Mathematics;
using UnityEngine;

[UpdateInGroup(typeof(FixedStepSimulationSystemGroup))]
[UpdateAfter(typeof(ClothSimulationSystem))]
[UpdateAfter(typeof(SoftBodySimulationSystem))]
[UpdateBefore(typeof(XPBDCollisionSystem))]
public partial struct CrossBodyCollisionSystem : ISystem
{
    private EntityQuery _clothQuery;
    private EntityQuery _softBodyQuery;

    // 迭代次数：每次迭代只处理"最深穿透"邻居，多迭代能逐个解开所有接触
    private const int NumIterations = 8;
    // 跨体摩擦：切向速度衰减比例
    private const float CrossBodyFriction = 0.2f;

    // 单次迭代最大修正比例：相对于 (rI+rJ)。
    // 1.0 表示一次最多推开一个"碰撞直径"距离；8 次迭代总上限 ~8*minDist 足以解开任何穿透
    private const float MaxCorrectionRatio = 1.0f;

    // 诊断：每 N 帧打一次日志，改成 0 可关闭
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

        var allClothEntities = _clothQuery.ToEntityArray(Allocator.Temp);
        var softBodyEntities = _softBodyQuery.ToEntityArray(Allocator.Temp);

        // 过滤掉 SkipCrossBodyCollision=true 的布料（例如升起/落下阶段的布料）
        // 这些布料本帧不参与跨体碰撞，避免正在被"锚点驱动"的布料和球之间产生怪异相互作用
        var clothEntitiesList = new NativeList<Entity>(allClothEntities.Length, Allocator.Temp);
        for (int i = 0; i < allClothEntities.Length; i++)
        {
            var cfg = state.EntityManager.GetComponentData<ClothSolverConfig>(allClothEntities[i]);
            if (!cfg.SkipCrossBodyCollision)
            {
                clothEntitiesList.Add(allClothEntities[i]);
            }
        }
        allClothEntities.Dispose();
        var clothEntities = clothEntitiesList.AsArray();

        int totalBodies = clothEntities.Length + softBodyEntities.Length;

        if (totalBodies < 2)
        {
            clothEntitiesList.Dispose();
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
            clothEntitiesList.Dispose();
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

        // === 2. 迭代：BuildHash → Resolve → Apply（全部并行）===
        // 并行块大小经验值：粒子数较少时用 32，能保持一定的局部性同时吃满多核
        const int kBatchSize = 32;
        JobHandle iterHandle = default;
        for (int iter = 0; iter < NumIterations; iter++)
        {
            // BuildHashJob 的 ParallelWriter 不支持 Clear，必须主线程先清空
            // 这里要等待上一轮的 iterHandle 完成后再 Clear（保证安全）
            iterHandle.Complete();
            hashMap.Clear();

            var buildHashJob = new BuildCrossBodyHashJob
            {
                GlobalPositions = globalPositions,
                HashMap = hashMap.AsParallelWriter(),
                CellSize = cellSize
            };
            iterHandle = buildHashJob.Schedule(totalParticles, kBatchSize);

            var resolveJob = new CrossBodyCollisionResolveJob
            {
                GlobalPositions = globalPositions,
                GlobalInvMasses = globalInvMasses,
                GlobalRadii = globalRadii,
                GlobalBodyIds = globalBodyIds,
                HashMap = hashMap,
                GlobalCorrections = globalCorrections,
                CellSize = cellSize,
                MaxCorrectionRatio = MaxCorrectionRatio
            };
            iterHandle = resolveJob.Schedule(totalParticles, kBatchSize, iterHandle);

            var applyJob = new ApplyGlobalCorrectionJob
            {
                GlobalPositions = globalPositions,
                GlobalCorrections = globalCorrections
            };
            iterHandle = applyJob.Schedule(totalParticles, kBatchSize, iterHandle);
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

        clothEntitiesList.Dispose();
        softBodyEntities.Dispose();
        offsets.Dispose();
        sizes.Dispose();
        radii.Dispose();
        frictions.Dispose();
    }

    // 彻底重写速度处理逻辑：
    // === 关键洞察 ===
    // 前面版本错在"手工维护速度"——先算 baseVel，再减去 velFromCorr，再加 relax*velFromCorr，
    // 再限幅。这一堆操作的语义混乱：当粒子被推开（corr 朝外），位置变了但 prev 没变，
    // 意味着 (newPos - prev)/dt 自然包含了"分离速度"，这才是物理上正确的速度。
    // 再去"velocityRelaxation * velFromCorr"反而会抵消掉这个分离趋势 → 看起来像吸引。
    //
    // === 新策略（标准 XPBD 做法）===
    // 1. 位置修正：newPos = oldPos + corr ⇒ 直接写 posBuf
    // 2. prev 不动（非摩擦情况下）⇒ (newPos - prev)/dt 天然产生分离速度
    // 3. 摩擦：把 prev 沿切向朝 newPos 拉近一个比例 friction
    //    ⇒ 结果是切向速度衰减 (1 - friction) 倍，法向速度完全保留
    // 4. velBuf 就写 (newPos - prev)/dt，不做任何限幅、松弛
    //
    // 这是最朴素、物理最正确的做法；不再有"粒子吸引"，穿透也能被正常顶开。
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

            // 1. 应用位置修正（posBuf[k] 可能已经被 PostSolve 写过 newPos，
            //    这里的 corr 是 "相对 originalPositions[gi] 的偏移"，要叠加到当前 posBuf 上）
            float3 oldPos = posBuf[k].Value;
            float3 newPos = oldPos + corr;
            posBuf[k] = new ParticlePosition { Value = newPos };

            // 2. 摩擦：prev 沿切向被朝 newPos 拉近（衰减切向速度）
            //    注意：法向上 prev 完全不动 → 法向速度 (newPos - prev)/dt 包含了分离趋势 → 不会"吸引"
            float3 prev = prevBuf[k].Value;
            if (friction > 0f)
            {
                float3 normal = corr / corrLen;
                float3 disp = newPos - prev;
                float tangDotN = math.dot(disp, normal);
                float3 tangent = disp - tangDotN * normal;
                float fric = math.saturate(friction);
                prev += tangent * fric;
                prevBuf[k] = new ParticlePrevPosition { Value = prev };
            }

            // 3. 速度 = (newPos - prev) / dt —— 纯物理推导，不做任何人为限幅/松弛
            //    这样：沿法向 corr → 产生朝外的分离速度 ✓
            //          沿切向摩擦 → 产生衰减的切向速度 ✓
            float3 finalVel = (newPos - prev) * invDt;
            velBuf[k] = new ParticleVelocity { Value = finalVel };
        }
    }
}
