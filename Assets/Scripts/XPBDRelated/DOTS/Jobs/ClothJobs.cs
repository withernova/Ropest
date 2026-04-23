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
    [ReadOnly] public NativeArray<ClothEdge> Edges;
    public NativeArray<float> Lambdas;
    public float Stiffness;
    public float Dt;

    public void Execute()
    {
        float alpha = Stiffness / (Dt * Dt);

        int edgeCount = Edges.Length;
        for (int i = 0; i < edgeCount; i++)
        {
            var edge = Edges[i];
            int id0 = edge.IndexA;
            int id1 = edge.IndexB;

            float w0 = InvMasses[id0];
            float w1 = InvMasses[id1];

            float3 diff = Positions[id0] - Positions[id1];
            float l = math.length(diff);

            if (l == 0f) continue;

            float3 gradC = diff / l;
            float l_rest = edge.RestLength;
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
/// 第二步：基于空间哈希的自碰撞检测与修正（IJob单线程，双向修正避免穿模）
/// 每个粒子只查询自身所在cell及相邻26个cell
/// 使用IJob保证双向修正的正确性，避免并行写入竞争
/// </summary>
[BurstCompile]
public struct ClothSelfCollisionJob : IJob
{
    public NativeArray<float3> Positions;
    [ReadOnly] public NativeArray<float> InvMasses;
    [ReadOnly] public NativeParallelMultiHashMap<int, int> HashMap;
    public float MinDistance;
    public float CellSize;
    public int Subdivision;     // 每行的顶点数（sub+1），用于判断相邻关系
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

// ============================================================
// Cloth 解析碰撞约束 Job（Burst编译，在SubStep内运行）
// 支持球体和Box碰撞体，纯数学运算，不依赖Unity Physics
// ============================================================
[BurstCompile]
public struct ClothAnalyticalCollisionJob : IJob
{
    public NativeArray<float3> Positions;
    [ReadOnly] public NativeArray<float> InvMasses;
    [ReadOnly] public NativeArray<AnalyticalColliderData> Colliders;
    public float ParticleRadius;    // 粒子碰撞半径
    public int NumParticles;

    public void Execute()
    {
        for (int i = 0; i < NumParticles; i++)
        {
            if (InvMasses[i] == 0f) continue;

            float3 pos = Positions[i];

            for (int c = 0; c < Colliders.Length; c++)
            {
                var col = Colliders[c];

                switch (col.Type)
                {
                    case AnalyticalColliderType.Sphere:
                        pos = ResolveSphereCollision(pos, col.Center, col.Radius, ParticleRadius);
                        break;
                    case AnalyticalColliderType.Box:
                        pos = ResolveBoxCollision(pos, col, ParticleRadius);
                        break;
                }
            }

            Positions[i] = pos;
        }
    }

    /// <summary>
    /// 球体碰撞：将粒子推到球面外
    /// </summary>
    static float3 ResolveSphereCollision(float3 particlePos, float3 sphereCenter, float sphereRadius, float particleRadius)
    {
        float3 diff = particlePos - sphereCenter;
        float dist = math.length(diff);
        float minDist = sphereRadius + particleRadius;

        if (dist < minDist && dist > 1e-8f)
        {
            float3 dir = diff / dist;
            particlePos = sphereCenter + dir * minDist;
        }
        else if (dist <= 1e-8f)
        {
            // 粒子在球心，随便推一个方向
            particlePos = sphereCenter + new float3(0, 1, 0) * minDist;
        }

        return particlePos;
    }

    /// <summary>
    /// Box碰撞：将粒子推到Box表面外
    /// 将粒子变换到Box局部空间，做AABB检测，再变换回世界空间
    /// </summary>
    static float3 ResolveBoxCollision(float3 particlePos, AnalyticalColliderData box, float particleRadius)
    {
        // 将粒子变换到Box局部空间
        float3 localPos = math.mul(box.InvRotation, particlePos - box.Center);

        // 扩展半尺寸（加上粒子半径）
        float3 expandedHalf = box.HalfExtents + particleRadius;

        // 检查是否在扩展Box内
        if (math.abs(localPos.x) < expandedHalf.x &&
            math.abs(localPos.y) < expandedHalf.y &&
            math.abs(localPos.z) < expandedHalf.z)
        {
            // 找到最近的面并推出
            float3 penetration = expandedHalf - math.abs(localPos);
            float3 sign = math.sign(localPos);

            // 选择穿透最浅的轴推出
            if (penetration.x <= penetration.y && penetration.x <= penetration.z)
            {
                localPos.x = sign.x * expandedHalf.x;
            }
            else if (penetration.y <= penetration.x && penetration.y <= penetration.z)
            {
                localPos.y = sign.y * expandedHalf.y;
            }
            else
            {
                localPos.z = sign.z * expandedHalf.z;
            }

            // 变换回世界空间
            particlePos = math.mul(box.Rotation, localPos) + box.Center;
        }

        return particlePos;
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
