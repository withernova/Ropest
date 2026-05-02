using Unity.Burst;
using Unity.Collections;
using Unity.Jobs;
using Unity.Mathematics;

// ============================================================================
// 布料专属 Job：构建空间哈希 + 自碰撞
// 说明：通用的 PreSolve / PostSolve / ResetLambda / DistanceConstraint /
//       AnalyticalCollision 已迁移至 XPBDCommonJobs.cs（与软体共用）。
// ============================================================================

[BurstCompile]
public struct BuildSpatialHashJob : IJob
{
    [ReadOnly] public NativeArray<float3> Positions;
    public NativeParallelMultiHashMap<int, int> HashMap;
    public float CellSize;
    public int NumParticles;

    public void Execute()
    {
        HashMap.Clear();
        float invCell = 1f / CellSize;
        for (int i = 0; i < NumParticles; i++)
        {
            int hash = SpatialHash(Positions[i], invCell);
            HashMap.Add(hash, i);
        }
    }

    static int SpatialHash(float3 pos, float invCell)
    {
        // 使用大质数做空间哈希，减少碰撞
        int x = (int)math.floor(pos.x * invCell);
        int y = (int)math.floor(pos.y * invCell);
        int z = (int)math.floor(pos.z * invCell);
        return x * 73856093 ^ y * 19349663 ^ z * 83492791;
    }
}

[BurstCompile]
public struct ClothSelfCollisionJob : IJob
{
    public NativeArray<float3> Positions;
    [ReadOnly] public NativeArray<float> InvMasses;
    [ReadOnly] public NativeParallelMultiHashMap<int, int> HashMap;
    public float MinDistance;
    public float CellSize;
    public int Subdivision;
    public int NumParticles;

    public void Execute()
    {
        float invCell = 1f / CellSize;
        float minDistSq = MinDistance * MinDistance;

        for (int i = 0; i < NumParticles; i++)
        {
            if (InvMasses[i] == 0f) continue;

            float3 posI = Positions[i];

            // 当前粒子所在的cell坐标
            int cx = (int)math.floor(posI.x * invCell);
            int cy = (int)math.floor(posI.y * invCell);
            int cz = (int)math.floor(posI.z * invCell);

            // 遍历3x3x3邻域（27个cell）
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
                            // 只处理 j > i 的对，避免重复处理
                            if (j <= i) continue;
                            if (InvMasses[j] == 0f) continue;

                            // 跳过网格拓扑上相邻的粒子
                            int diff = math.abs(j - i);
                            if (diff == 1 || diff == Subdivision || diff == Subdivision - 1 || diff == Subdivision + 1)
                                continue;

                            float3 diff3 = Positions[i] - Positions[j];
                            float distSq = math.lengthsq(diff3);

                            if (distSq < minDistSq && distSq > 1e-12f)
                            {
                                float dist = math.sqrt(distSq);
                                float3 dir = diff3 / dist;
                                float overlap = MinDistance - dist;

                                float w0 = InvMasses[i];
                                float w1 = InvMasses[j];
                                float wSum = w0 + w1;
                                if (wSum < 1e-8f) continue;

                                // 双向修正：两个粒子都被推开
                                float3 corrI = (w0 / wSum) * overlap * dir;
                                float3 corrJ = -(w1 / wSum) * overlap * dir;

                                Positions[i] = Positions[i] + corrI;
                                Positions[j] = Positions[j] + corrJ;
                            }
                        } while (HashMap.TryGetNextValue(out j, ref it));
                    }
                }
            }
        }
    }
}
