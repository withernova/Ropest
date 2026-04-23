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
//   迭代 N 次：BuildCrossBodyHashJob → CrossBodyCollisionResolveJob → ApplyGlobalCorrectionJob →
//   主线程 Scatter（直接循环写回 Buffer）
//
// 注意：Gather 和 Scatter 改为在主线程完成（避免 ECS 对同一 ComponentType 的
//      Buffer 在 Schedule 期间多次访问引发的 AtomicSafetyHandle 冲突）。
//      因此这里只保留 Burst 化的三个核心计算 Job。
// ============================================================

/// <summary>
/// 构建全局空间哈希表（所有实体的粒子共用一张表）
/// </summary>
[BurstCompile]
public struct BuildCrossBodyHashJob : IJob
{
    [ReadOnly] public NativeArray<float3> GlobalPositions;
    public NativeParallelMultiHashMap<int, int> HashMap;
    public float CellSize;
    public int NumParticles;

    public void Execute()
    {
        HashMap.Clear();
        float invCell = 1f / CellSize;
        for (int i = 0; i < NumParticles; i++)
        {
            int hash = SpatialHash(GlobalPositions[i], invCell);
            HashMap.Add(hash, i);
        }
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
/// 跨体碰撞解算：只处理 BodyId 不同的粒子对，双向推开，按质量反比分配位移
/// 使用 IJob 单线程执行，保证双向修正的正确性（避免并行写竞争）
/// 修正结果写到 GlobalCorrections（进入时会先自清零）
/// </summary>
[BurstCompile]
public struct CrossBodyCollisionResolveJob : IJob
{
    [ReadOnly] public NativeArray<float3> GlobalPositions;
    [ReadOnly] public NativeArray<float> GlobalInvMasses;
    [ReadOnly] public NativeArray<float> GlobalRadii;
    [ReadOnly] public NativeArray<int> GlobalBodyIds;
    [ReadOnly] public NativeParallelMultiHashMap<int, int> HashMap;

    public NativeArray<float3> GlobalCorrections;

    public float CellSize;
    public int NumParticles;

    public void Execute()
    {
        // 进入时清零修正量
        for (int k = 0; k < NumParticles; k++)
            GlobalCorrections[k] = float3.zero;

        float invCell = 1f / CellSize;

        for (int i = 0; i < NumParticles; i++)
        {
            float wI = GlobalInvMasses[i];
            float3 posI = GlobalPositions[i];
            int bodyI = GlobalBodyIds[i];
            float rI = GlobalRadii[i];

            int cx = (int)math.floor(posI.x * invCell);
            int cy = (int)math.floor(posI.y * invCell);
            int cz = (int)math.floor(posI.z * invCell);

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
                            if (j <= i) continue;
                            // 只处理不同 Body 之间的碰撞（同 Body 自碰撞由其自身系统处理）
                            if (GlobalBodyIds[j] == bodyI) continue;

                            float wJ = GlobalInvMasses[j];
                            if (wI + wJ < 1e-8f) continue; // 两个都是固定点，无法处理

                            // === 跨体修正权重分配 ===
                            // 不使用"真实 invMass 反比"——因为 Cloth 质量基于面积、SoftBody 基于体积，
                            // 两套量纲完全不同（实测 Cloth invMass≈0.89，SoftBody invMass≈200），
                            // 会导致"轻的一方吸收 99% 修正，重的一方几乎不动"的严重失衡。
                            // 改为"是否可动"的二值分配：双方都能动 → 50/50 均分；
                            // 一方固定 → 另一方承担全部修正。
                            float shareI, shareJ;
                            if (wI > 0f && wJ > 0f)
                            {
                                shareI = 0.5f;
                                shareJ = 0.5f;
                            }
                            else if (wI > 0f)
                            {
                                shareI = 1.0f;
                                shareJ = 0f;
                            }
                            else
                            {
                                shareI = 0f;
                                shareJ = 1.0f;
                            }

                            float rJ = GlobalRadii[j];
                            float minDist = rI + rJ;
                            float minDistSq = minDist * minDist;

                            float3 diff = GlobalPositions[i] - GlobalPositions[j];
                            float distSq = math.lengthsq(diff);
                            if (distSq >= minDistSq) continue;

                            if (distSq < 1e-12f)
                            {
                                // 完全重合时随便选个方向
                                float3 fallback = new float3(0, 1, 0);
                                GlobalCorrections[i] += shareI * minDist * fallback;
                                GlobalCorrections[j] -= shareJ * minDist * fallback;
                                continue;
                            }

                            float dist = math.sqrt(distSq);
                            float3 dir = diff / dist;
                            float overlap = minDist - dist;

                            GlobalCorrections[i] += shareI * overlap * dir;
                            GlobalCorrections[j] -= shareJ * overlap * dir;

                        } while (HashMap.TryGetNextValue(out j, ref it));
                    }
                }
            }
        }
    }
}
