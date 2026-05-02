using Unity.Burst;
using Unity.Collections;
using Unity.Jobs;
using Unity.Mathematics;


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
    public float MaxCorrectionRatio;

    public void Execute(int i)
    {
        float wI = GlobalInvMasses[i];
        if (wI <= 0f)
        {
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

        // === 重写：单次迭代只处理"最深穿透"的那个邻居 ===
        // 为什么：之前尝试过"每轴取最大绝对值"→ 方向混乱导致发散；"加权平均法线"→ 
        //        多邻居方向冲突时幅度不足导致穿透/吸引。都不如最朴素的做法：
        // 本次迭代只采纳单个最深穿透邻居的推开向量，多邻居通过多次迭代自然收敛。
        // 这是 Unity Physics / PhysX / Bullet 都使用的"sequential impulse"思路的位置版。
        // —— 单邻居方向是真实接触法向，没有任何"合成误差"；
        // —— 多邻居冲突？下一次迭代会选择到当前最深的那个；
        // —— 配合 NumIterations=6+，足以让所有接触都被逐个解开。

        float3 deepestCorr = float3.zero;
        float deepestOverlap = 0f;
        float deepestMinDist = 0f;

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
                            // 退化：两点重合。用稳定的确定性方向避免同向弹飞
                            int seed = (i - j);
                            dir = math.normalize(new float3(
                                ((seed * 12.9898f) % 1f) - 0.5f,
                                ((seed * 78.2330f) % 1f) - 0.5f + 0.1f,
                                ((seed * 37.7190f) % 1f) - 0.5f));
                            overlap = minDist;
                        }
                        else
                        {
                            float dist = math.sqrt(distSq);
                            dir = diff / dist;
                            overlap = minDist - dist;
                        }

                        // 只记录最深的那个
                        if (overlap > deepestOverlap)
                        {
                            deepestOverlap = overlap;
                            deepestCorr = dir * (shareI * overlap);
                            deepestMinDist = minDist;
                        }

                    } while (HashMap.TryGetNextValue(out j, ref it));
                }
            }
        }

        if (deepestOverlap <= 0f)
        {
            GlobalCorrections[i] = float3.zero;
            return;
        }

        // 单次迭代幅度限幅（防止初始大穿透一步跳太远）
        float maxLen = MaxCorrectionRatio * deepestMinDist;
        float corrLen = math.length(deepestCorr);
        if (corrLen > maxLen && corrLen > 1e-8f)
        {
            deepestCorr *= (maxLen / corrLen);
        }

        GlobalCorrections[i] = deepestCorr;
    }
}

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
