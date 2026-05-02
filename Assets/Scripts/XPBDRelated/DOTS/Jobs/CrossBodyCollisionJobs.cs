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
