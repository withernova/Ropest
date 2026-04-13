using Unity.Burst;
using Unity.Collections;
using Unity.Jobs;
using Unity.Mathematics;

// ============================================================
// Cloth PreSolve Job - 重力积分（可并行）
// ============================================================
[BurstCompile]
public struct ClothPreSolveJob : IJobParallelFor
{
    public NativeArray<float3> Positions;
    public NativeArray<float3> PrevPositions;
    public NativeArray<float3> Velocities;
    [ReadOnly] public NativeArray<float> InvMasses;
    public float3 Gravity;
    public float Dt;

    public void Execute(int i)
    {
        if (InvMasses[i] == 0f) return;

        float3 vel = Velocities[i];
        vel += Gravity * Dt;

        PrevPositions[i] = Positions[i];
        Velocities[i] = vel;
        Positions[i] += vel * Dt;
    }
}

// ============================================================
// Cloth Lambda重置 Job
// ============================================================
[BurstCompile]
public struct ClothResetLambdaJob : IJobParallelFor
{
    public NativeArray<float> Lambdas;

    public void Execute(int i)
    {
        Lambdas[i] = 0f;
    }
}

// ============================================================
// Cloth DistanceConstraint Job（顺序依赖，使用IJob + Burst）
// 注意：如果需要更高并行度，可以用图着色分组后IJobParallelFor
// ============================================================
[BurstCompile]
public struct ClothDistanceConstraintJob : IJob
{
    public NativeArray<float3> Positions;
    [ReadOnly] public NativeArray<float> InvMasses;
    [ReadOnly] public NativeArray<int> EdgeIndexA;
    [ReadOnly] public NativeArray<int> EdgeIndexB;
    [ReadOnly] public NativeArray<float> RestLengths;
    public NativeArray<float> Lambdas;
    public float Stiffness;
    public float Dt;

    public void Execute()
    {
        float alpha = Stiffness / (Dt * Dt);

        for (int i = 0; i < EdgeIndexA.Length; i++)
        {
            int id0 = EdgeIndexA[i];
            int id1 = EdgeIndexB[i];

            float w0 = InvMasses[id0];
            float w1 = InvMasses[id1];

            float3 diff = Positions[id0] - Positions[id1];
            float l = math.length(diff);

            if (l == 0f) continue;

            float3 gradC = diff / l;
            float l_rest = RestLengths[i];
            float C = l - l_rest;
            float wTot = (w1 + w0) * math.lengthsq(gradC);

            float deltaLambda = -(C + Lambdas[i] * alpha) / (wTot + alpha);
            Lambdas[i] += deltaLambda;

            Positions[id0] += deltaLambda * w0 * gradC;
            Positions[id1] -= deltaLambda * w1 * gradC;
        }
    }
}

// ============================================================
// Cloth 自碰撞约束 - 空间哈希加速版（Burst编译）
// 将空间划分为网格，每个粒子只检测相邻cell中的粒子
// 复杂度从O(n²)降到接近O(n)
// ============================================================

/// <summary>
/// 第一步：构建空间哈希表（IJob，单线程写入HashMap）
/// </summary>
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

/// <summary>
/// 第二步：基于空间哈希的自碰撞检测与修正（IJobParallelFor，并行处理每个粒子）
/// 每个粒子只查询自身所在cell及相邻26个cell
/// </summary>
[BurstCompile]
public struct ClothSelfCollisionJob : IJobParallelFor
{
    [NativeDisableParallelForRestriction]
    public NativeArray<float3> Positions;
    [ReadOnly] public NativeArray<float> InvMasses;
    [ReadOnly] public NativeParallelMultiHashMap<int, int> HashMap;
    public float MinDistance;
    public float CellSize;
    public int Subdivision;     // 每行的顶点数（sub+1），用于判断相邻关系

    public void Execute(int i)
    {
        if (InvMasses[i] == 0f) return;

        float3 posI = Positions[i];
        float invCell = 1f / CellSize;
        float minDistSq = MinDistance * MinDistance;

        // 当前粒子所在的cell坐标
        int cx = (int)math.floor(posI.x * invCell);
        int cy = (int)math.floor(posI.y * invCell);
        int cz = (int)math.floor(posI.z * invCell);

        float3 totalCorrection = float3.zero;
        int correctionCount = 0;

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

                        float3 diff3 = posI - Positions[j];
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

                            // 累积修正量（并行安全：只修改自己的位置）
                            totalCorrection += (w0 / wSum) * overlap * dir;
                            correctionCount++;
                        }
                    } while (HashMap.TryGetNextValue(out j, ref it));
                }
            }
        }

        // 应用累积修正
        if (correctionCount > 0)
        {
            Positions[i] = posI + totalCorrection;
        }
    }
}

// ============================================================
// Cloth PostSolve Job - 速度更新 + 阻尼（可并行）
// ============================================================
[BurstCompile]
public struct ClothPostSolveJob : IJobParallelFor
{
    public NativeArray<float3> Positions;
    [ReadOnly] public NativeArray<float3> PrevPositions;
    public NativeArray<float3> Velocities;
    [ReadOnly] public NativeArray<float> InvMasses;
    public float OneOverDt;
    public float Dt;
    public float Damping; // 阻尼系数

    public void Execute(int i)
    {
        if (InvMasses[i] == 0f) return;

        // 标准速度阻尼：直接衰减位移量（等价于 vel *= (1 - damping)）
        // Damping=0 无阻尼，Damping=1 完全静止
        // 推荐范围 0.01~0.05
        if (Damping > 0f)
        {
            float3 displacement = Positions[i] - PrevPositions[i];
            // 将位移缩小，等价于速度衰减
            Positions[i] = PrevPositions[i] + displacement * (1f - Damping);
        }

        Velocities[i] = (Positions[i] - PrevPositions[i]) * OneOverDt;
    }
}
