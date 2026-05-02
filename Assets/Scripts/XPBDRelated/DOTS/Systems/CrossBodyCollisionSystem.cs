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

        // === 分配全局扇平数组 ===
        var globalPositions = new NativeArray<float3>(totalParticles, Allocator.TempJob, NativeArrayOptions.UninitializedMemory);
        var globalPrevPositions = new NativeArray<float3>(totalParticles, Allocator.TempJob, NativeArrayOptions.UninitializedMemory);
        var globalInvMasses = new NativeArray<float>(totalParticles, Allocator.TempJob, NativeArrayOptions.UninitializedMemory);
        var globalRadii = new NativeArray<float>(totalParticles, Allocator.TempJob, NativeArrayOptions.UninitializedMemory);
        var globalBodyIds = new NativeArray<int>(totalParticles, Allocator.TempJob, NativeArrayOptions.UninitializedMemory);
        // 是否为表面粒子（1=表面，参与跨体碰撞；0=内部，跳过）
        var globalIsSurface = new NativeArray<byte>(totalParticles, Allocator.TempJob, NativeArrayOptions.UninitializedMemory);
        var globalCorrections = new NativeArray<float3>(totalParticles, Allocator.TempJob, NativeArrayOptions.ClearMemory);

        // === 每个 body 的粒子质心（每轮迭代重新计算）===
        // 质心用于"远离对方 body 整体中心"的方向参考：当一个粒子严重穿进对方球内时，
        // 基于粒子-粒子连线的 dir 可能指向对方更深处（错方向），但基于粒子-对方质心的 dir
        // 永远指向球外，这才是稳健的分离方向。
        var bodyCenters = new NativeArray<float3>(totalBodies, Allocator.TempJob, NativeArrayOptions.UninitializedMemory);
        // 每个 body 的"有效外接半径"（表面粒子到质心的最大距离），用于 body 级 shell 碰撞约束
        var bodyRadii = new NativeArray<float>(totalBodies, Allocator.TempJob, NativeArrayOptions.UninitializedMemory);

        // === 1. 主线程收集：把所有实体的 Buffer 数据拷到全局数组 ===
        for (int i = 0; i < numCloth; i++)
        {
            var entity = clothEntities[i];
            var posBuf = state.EntityManager.GetBuffer<ParticlePosition>(entity, true);
            var prevBuf = state.EntityManager.GetBuffer<ParticlePrevPosition>(entity, true);
            var invMassBuf = state.EntityManager.GetBuffer<ParticleInvMass>(entity, true);
            // 布料全部粒子视为表面（若没有 ParticleSurfaceFlag Buffer 也按表面处理）
            bool hasFlagBuf = state.EntityManager.HasBuffer<ParticleSurfaceFlag>(entity);
            DynamicBuffer<ParticleSurfaceFlag> flagBuf = default;
            if (hasFlagBuf) flagBuf = state.EntityManager.GetBuffer<ParticleSurfaceFlag>(entity, true);
            int num = sizes[i];
            int off = offsets[i];
            float r = radii[i];
            int bid = i;
            for (int k = 0; k < num; k++)
            {
                globalPositions[off + k] = posBuf[k].Value;
                globalPrevPositions[off + k] = prevBuf[k].Value;
                globalInvMasses[off + k] = invMassBuf[k].Value;
                globalRadii[off + k] = r;
                globalBodyIds[off + k] = bid;
                // 布料默认全部是表面粒子
                globalIsSurface[off + k] = hasFlagBuf ? flagBuf[k].Value : (byte)1;
            }
        }
        for (int i = 0; i < numSoft; i++)
        {
            int idx = numCloth + i;
            var entity = softBodyEntities[i];
            var posBuf = state.EntityManager.GetBuffer<ParticlePosition>(entity, true);
            var prevBuf = state.EntityManager.GetBuffer<ParticlePrevPosition>(entity, true);
            var invMassBuf = state.EntityManager.GetBuffer<ParticleInvMass>(entity, true);
            // 软体：只有表面粒子（ParticleSurfaceFlag.Value==1）参与跨体碰撞
            bool hasFlagBuf = state.EntityManager.HasBuffer<ParticleSurfaceFlag>(entity);
            DynamicBuffer<ParticleSurfaceFlag> flagBuf = default;
            if (hasFlagBuf) flagBuf = state.EntityManager.GetBuffer<ParticleSurfaceFlag>(entity, true);
            int num = sizes[idx];
            int off = offsets[idx];
            float r = radii[idx];
            int bid = idx;
            for (int k = 0; k < num; k++)
            {
                globalPositions[off + k] = posBuf[k].Value;
                globalPrevPositions[off + k] = prevBuf[k].Value;
                globalInvMasses[off + k] = invMassBuf[k].Value;
                globalRadii[off + k] = r;
                globalBodyIds[off + k] = bid;
                // 没有标记 Buffer 时退化为全表面（向后兼容，但建议旧软体重建一次）
                globalIsSurface[off + k] = hasFlagBuf ? flagBuf[k].Value : (byte)1;
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

            // 每轮都重算 body 质心（迭代过程中粒子在动，质心也随之变化）
            var centerJob = new ComputeBodyCentersJob
            {
                GlobalPositions = globalPositions,
                GlobalBodyIds = globalBodyIds,
                GlobalIsSurface = globalIsSurface,
                TotalParticles = totalParticles,
                BodyCenters = bodyCenters,
                BodyRadii = bodyRadii
            };
            iterHandle = centerJob.Schedule();

            var buildHashJob = new BuildCrossBodyHashJob
            {
                GlobalPositions = globalPositions,
                GlobalIsSurface = globalIsSurface,
                HashMap = hashMap.AsParallelWriter(),
                CellSize = cellSize
            };
            iterHandle = buildHashJob.Schedule(totalParticles, kBatchSize, iterHandle);

            var resolveJob = new CrossBodyCollisionResolveJob
            {
                GlobalPositions = globalPositions,
                GlobalPrevPositions = globalPrevPositions,
                GlobalInvMasses = globalInvMasses,
                GlobalRadii = globalRadii,
                GlobalBodyIds = globalBodyIds,
                GlobalIsSurface = globalIsSurface,
                BodyCenters = bodyCenters,
                BodyRadii = bodyRadii,
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
        globalPrevPositions.Dispose();
        originalPositions.Dispose();
        globalInvMasses.Dispose();
        globalRadii.Dispose();
        globalBodyIds.Dispose();
        globalIsSurface.Dispose();
        globalCorrections.Dispose();
        hashMap.Dispose();
        bodyCenters.Dispose();
        bodyRadii.Dispose();

        clothEntitiesList.Dispose();
        softBodyEntities.Dispose();
        offsets.Dispose();
        sizes.Dispose();
        radii.Dispose();
        frictions.Dispose();
    }

    // === 方案 A：显式拆分法向/切向并钳制"朝向对方"的入射速度 ===
    //
    // 旧做法 finalVel = (newPos - prev)/dt 的问题：
    //   prev 是帧初位置，当入射速度巨大时 (oldPos - prev) 方向朝对方，其分量可能远大于 corr，
    //   导致 finalVel 仍朝对方 → 下一帧 PreSolve 又穿进去 → 视觉上"相互吸引/黏着"。
    //
    // 新做法（基于 PostSolve 已写好的 velBuf）：
    //   1. 取出 PostSolve 的速度 oldVel，沿 normal=corr/|corr| 分解为 vn（法向）、vt（切向）
    //   2. 法向：如果粒子仍朝对方（vn < 0，朝着 -normal 方向=对方方向），钳制到至少分离速度 corrLen/dt；
    //      如果已经在远离（vn >= 0），保留原速度（不要"加"分离速度，避免加速抛飞）
    //   3. 切向：按摩擦系数衰减 vt *= (1 - friction)
    //   4. finalVel = vt + vnNew * normal
    //   5. 反推 prev = newPos - finalVel * dt，让 (newPos - prev)/dt 与 finalVel 自洽，
    //      下一帧 PreSolve 读取 vel 后位置积分正确，不会再撞回去。
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

            // 1. 应用位置修正
            float3 oldPos = posBuf[k].Value;
            float3 newPos = oldPos + corr;
            posBuf[k] = new ParticlePosition { Value = newPos };

            // 2. 以 PostSolve 写好的速度为基准做法向/切向重组
            float3 normal = corr / corrLen;   // 指向远离对方（分离方向）
            float3 oldVel = velBuf[k].Value;
            float vn = math.dot(oldVel, normal);           // 法向分量：>0 分离，<0 朝对方
            float3 vt = oldVel - vn * normal;              // 切向分量

            // 3. 法向钳制：
            //    - 如果 vn < 0（朝对方）→ 抹掉入射速度，替换为 corrLen/dt 的分离速度
            //    - 如果 0 <= vn < corrLen/dt → 提升到 corrLen/dt 保证能分离
            //    - 如果 vn 已经 >= corrLen/dt → 保留原值
            float vnSeparation = corrLen * invDt;
            float vnNew = math.max(vn, vnSeparation);

            // 4. 切向摩擦：标准 (1 - friction) 衰减
            float fric = math.saturate(friction);
            vt *= (1f - fric);

            float3 finalVel = vt + vnNew * normal;

            // 5. 反推 prev 保证 (newPos - prev)/dt == finalVel，
            //    这样下一帧 PreSolve 写入的 PrevPositions = Positions 之前，
            //    Damping 等对 prev 的假设依然成立，且不存在"残留的入射速度"。
            float3 newPrev = newPos - finalVel * dt;
            prevBuf[k] = new ParticlePrevPosition { Value = newPrev };
            velBuf[k] = new ParticleVelocity { Value = finalVel };
        }
    }
}
