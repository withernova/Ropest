using Unity.Burst;
using Unity.Collections;
using Unity.Jobs;
using Unity.Mathematics;

// ============================================================
// 跨体碰撞（Cross-Body Collision）所需的 Job 集合
// 目标：让不同实体（Cloth、SoftBody）的粒子之间互相推开，
// 形成真实的双向形变耦合（例如软体球砸到布料上，球被顶出凹陷，布料也被顶起）
//
// 工作流程（与 CrossBodyCollisionSystem 配合）：
//   主线程 Gather（直接循环拷贝 Buffer）→
//   迭代 N 次：[主线程 HashMap.Clear] → BuildCrossBodyHashJob(并行) →
//                CrossBodyCollisionResolveJob(并行) → ApplyGlobalCorrectionJob(并行) →
//   主线程 Scatter（直接循环写回 Buffer）
//
// 并行化策略（2026-04 优化）：
//   - BuildCrossBodyHashJob: IJobParallelFor + NativeParallelMultiHashMap.ParallelWriter
//   - CrossBodyCollisionResolveJob: IJobParallelFor，每个粒子 i 只写 GlobalCorrections[i]
//     （不再写 j），通过不跳过 j<=i 让 (i,j) 和 (j,i) 两对独立计算，
//     利用 diff 方向对称性得到与原版等价的双向修正结果
//   - ApplyGlobalCorrectionJob: IJobParallelFor
//
// 注意：Gather 和 Scatter 仍然在主线程完成（避免 ECS 对同一 ComponentType 的
//      Buffer 在 Schedule 期间多次访问引发的 AtomicSafetyHandle 冲突）。
// ============================================================

/// <summary>
/// 构建全局空间哈希表（所有实体的粒子共用一张表）
/// 并行版：每个线程只负责一部分粒子的插入，用 ParallelWriter 保证线程安全
/// 注意：调用前必须由主线程先 HashMap.Clear()（ParallelWriter 不支持 Clear）
/// </summary>
[BurstCompile]
public struct BuildCrossBodyHashJob : IJobParallelFor
{
    [ReadOnly] public NativeArray<float3> GlobalPositions;
    public NativeParallelMultiHashMap<int, int>.ParallelWriter HashMap;
    public float CellSize;

    public void Execute(int i)
    {
        float invCell = 1f / CellSize;
        int hash = SpatialHash(GlobalPositions[i], invCell);
        HashMap.Add(hash, i);
    }

    public static int SpatialHash(float3 pos, float invCell)
    {
        int x = (int)math.floor(pos.x * invCell);
        int y = (int)math.floor(pos.y * invCell);
        int z = (int)math.floor(pos.z * invCell);
        return x * 73856093 ^ y * 19349663 ^ z * 83492791;
    }
}

/// <summary>
/// 跨体碰撞解算 - 并行版（已改进为"加权平均 + 夹紧"，避免局部粒子被多邻居累加弹飞）
/// 每个粒子 i 在自己的线程里遍历 27 个邻居 cell，
/// 用"按穿透深度加权平均"替代原本的"线性累加"，得到一个更接近真实"被一群粒子顶住"的
/// 平均位移，避免修正量被 N 个重叠邻居放大 N 倍，从而消除局部突出/抽搐。
///
/// 数学：
///   accumCorr = Σ (overlap_k * overlap_k * dir_k)   // 用 overlap^2 做权重
///   accumW    = Σ overlap_k
///   finalCorr = accumCorr / accumW
/// 直观上等价于"取各邻居推开方向的穿透加权重心"，对一群穿透深度相近的邻居，
/// 结果 ≈ 一个邻居单独推的距离，避免了"多邻居 → 总位移爆表"的抽搐现象。
///
/// 此外修正还做了安全夹紧：|finalCorr| ≤ MaxCorrectionRatio * (rI + rJ_max)，
/// 防止单帧推距过大引起视觉上的"跳跃"。
///
/// 关键属性：
///   - NativeDisableParallelForRestriction 不需要（每线程只写 index=i 的 correction）
///   - 仅 GlobalCorrections 为 WriteOnly-per-index，符合 IJobParallelFor 的安全规则
/// </summary>
[BurstCompile]
public struct CrossBodyCollisionResolveJob : IJobParallelFor
{
    [ReadOnly] public NativeArray<float3> GlobalPositions;
    [ReadOnly] public NativeArray<float> GlobalInvMasses;
    [ReadOnly] public NativeArray<float> GlobalRadii;
    [ReadOnly] public NativeArray<int> GlobalBodyIds;
    [ReadOnly] public NativeParallelMultiHashMap<int, int> HashMap;

    public NativeArray<float3> GlobalCorrections;

    public float CellSize;

    // 单次迭代单个粒子的最大修正幅度 = MaxCorrectionRatio * (rI + rJ)
    // 1.0 表示最多一次性推开 1 倍"碰撞直径"的距离，足够把穿透解开但不会跳跃
    public float MaxCorrectionRatio;

    public void Execute(int i)
    {
        float wI = GlobalInvMasses[i];
        if (wI <= 0f)
        {
            // 固定粒子：不累积任何修正
            GlobalCorrections[i] = float3.zero;
            return;
        }

        float3 posI = GlobalPositions[i];
        int bodyI = GlobalBodyIds[i];
        float rI = GlobalRadii[i];

        float invCell = 1f / CellSize;
        int cx = (int)math.floor(posI.x * invCell);
        int cy = (int)math.floor(posI.y * invCell);
        int cz = (int)math.floor(posI.z * invCell);

        // 加权平均累加器
        float3 weightedCorr = float3.zero;
        float totalW = 0f;
        float maxMinDist = 0f; // 用于最终 clamp

        for (int dx = -1; dx <= 1; dx++)
        {
            for (int dy = -1; dy <= 1; dy++)
            {
                for (int dz = -1; dz <= 1; dz++)
                {
                    int hash = (cx + dx) * 73856093 ^ (cy + dy) * 19349663 ^ (cz + dz) * 83492791;

                    if (!HashMap.TryGetFirstValue(hash, out int j, out var it)) continue;

                    do
                    {
                        if (j == i) continue;
                        if (GlobalBodyIds[j] == bodyI) continue;

                        float wJ = GlobalInvMasses[j];
                        if (wI + wJ < 1e-8f) continue;

                        // 跨体修正权重：双方都能动 → 本侧 0.5；对方固定 → 本侧 1.0
                        float shareI = (wJ > 0f) ? 0.5f : 1.0f;

                        float rJ = GlobalRadii[j];
                        float minDist = rI + rJ;
                        float minDistSq = minDist * minDist;

                        float3 diff = posI - GlobalPositions[j];
                        float distSq = math.lengthsq(diff);
                        if (distSq >= minDistSq) continue;

                        float3 dir;
                        float overlap;

                        if (distSq < 1e-12f)
                        {
                            // 完全重合 fallback：使用稳定的 +Y 方向，不再用 (i-j) 这种
                            // 依赖粒子编号的方向，避免横向抽动。
                            // i, j 同时命中 (j,i) 时方向也都是 +Y，不会互相抵消吗？
                            // → 不会：(j,i) 那一边 diff = posJ - posI ≈ 0，它也会走这条路径，
                            //   但它只写自己（j）的 correction，不会影响 i 的方向选择。
                            //   实际仅在完全重合的极端情况才会进入此分支，可忽略方向一致性。
                            dir = new float3(0f, 1f, 0f);
                            overlap = minDist;
                        }
                        else
                        {
                            float dist = math.sqrt(distSq);
                            dir = diff / dist;
                            overlap = minDist - dist;
                        }

                        // 用 overlap 作权重：穿透浅的邻居贡献小，深的贡献大
                        // 这样密集邻居（浅穿透）被几何平均后不会爆表
                        float w = overlap;
                        weightedCorr += shareI * overlap * dir * w;
                        totalW += w;

                        if (minDist > maxMinDist) maxMinDist = minDist;

                    } while (HashMap.TryGetNextValue(out j, ref it));
                }
            }
        }

        if (totalW < 1e-8f)
        {
            GlobalCorrections[i] = float3.zero;
            return;
        }

        // 加权平均：得到一个等效的"合力方向 × 平均穿透深度"的修正量
        float3 finalCorr = weightedCorr / totalW;

        // 夹紧到 MaxCorrectionRatio * minDist，防止极端情况（如初始穿透太深）一帧跳太远
        float maxLen = MaxCorrectionRatio * maxMinDist;
        float corrLen = math.length(finalCorr);
        if (corrLen > maxLen && corrLen > 1e-8f)
        {
            finalCorr *= (maxLen / corrLen);
        }

        GlobalCorrections[i] = finalCorr;
    }
}

/// <summary>
/// 把本次迭代累积的 Corrections 应用到 GlobalPositions - 并行版
/// 原版串行循环改为 IJobParallelFor，每索引独立写回
/// </summary>
[BurstCompile]
public struct ApplyGlobalCorrectionJob : IJobParallelFor
{
    public NativeArray<float3> GlobalPositions;
    public NativeArray<float3> GlobalCorrections;

    public void Execute(int i)
    {
        GlobalPositions[i] += GlobalCorrections[i];
        GlobalCorrections[i] = float3.zero;
    }
}
