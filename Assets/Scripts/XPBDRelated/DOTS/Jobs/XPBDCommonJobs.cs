using Unity.Burst;
using Unity.Collections;
using Unity.Jobs;
using Unity.Mathematics;

// ============================================================================
// XPBD 通用流水线 Job（布料与软体共用）
// 说明：布料与软体的 PreSolve / PostSolve / ResetLambda / DistanceConstraint /
//       AnalyticalCollision 原本逻辑几乎完全一致（仅 Edge/Lambda 结构体不同），
//       现通过统一 XPBDEdge / XPBDDistanceLambda 后合并成本文件这一组共享 Job。
//       专属逻辑（布料自碰撞、软体体积约束）仍放在各自的 ClothJobs / SoftBodyJobs。
// ============================================================================

/// <summary>
/// 外力积分与预测位置（重力 * dt 累加到速度，位置积分到 Predicted）。
/// 对应旧版：ClothPreSolveJob / SoftBodyPreSolveJob。
/// </summary>
[BurstCompile]
public struct XPBDPreSolveJob : IJobParallelFor
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

/// <summary>
/// 将一个 float Buffer（Lambda 累积量）重置为 0。
/// 对应旧版：ClothResetLambdaJob / SoftBodyResetDistanceLambdaJob / SoftBodyResetVolumeLambdaJob。
/// </summary>
[BurstCompile]
public struct XPBDResetFloatBufferJob : IJobParallelFor
{
    public NativeArray<float> Values;

    public void Execute(int i)
    {
        Values[i] = 0f;
    }
}

/// <summary>
/// XPBD 距离约束（所有 Edge 的等距离约束），串行 IJob 保证同一粒子写入无竞争。
/// 对应旧版：ClothDistanceConstraintJob / SoftBodyDistanceConstraintJob。
/// </summary>
[BurstCompile]
public struct XPBDDistanceConstraintJob : IJob
{
    public NativeArray<float3> Positions;
    [ReadOnly] public NativeArray<float> InvMasses;
    [ReadOnly] public NativeArray<XPBDEdge> Edges;
    public NativeArray<float> Lambdas;
    public float Stiffness; // compliance
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

/// <summary>
/// XPBD 距离约束（图着色并行版本）。
/// 调用方必须保证 [ColorStart, ColorStart+ColorCount) 范围内的边两两不共享粒子，
/// 否则会出现写竞争。每个颜色组调度一次 IJobParallelFor，色间串行、色内并行。
///
/// NativeDisableParallelForRestriction：允许在 ParallelFor 里按任意索引写 Positions，
/// 因为同色组内各边写入的粒子下标互不相同，由图着色算法保证无冲突。
/// </summary>
[BurstCompile]
public struct XPBDDistanceConstraintColoredJob : IJobParallelFor
{
    [NativeDisableParallelForRestriction] public NativeArray<float3> Positions;
    [ReadOnly] public NativeArray<float> InvMasses;
    [ReadOnly] public NativeArray<XPBDEdge> Edges;
    [NativeDisableParallelForRestriction] public NativeArray<float> Lambdas;
    public float Stiffness;
    public float Dt;
    public int ColorStart;

    public void Execute(int localIndex)
    {
        int i = ColorStart + localIndex;
        var edge = Edges[i];
        int id0 = edge.IndexA;
        int id1 = edge.IndexB;

        float w0 = InvMasses[id0];
        float w1 = InvMasses[id1];

        float3 diff = Positions[id0] - Positions[id1];
        float l = math.length(diff);
        if (l == 0f) return;

        float3 gradC = diff / l;
        float l_rest = edge.RestLength;
        float C = l - l_rest;
        float wTot = (w1 + w0) * math.lengthsq(gradC);

        float alpha = Stiffness / (Dt * Dt);
        float deltaLambda = -(C + Lambdas[i] * alpha) / (wTot + alpha);
        Lambdas[i] += deltaLambda;

        Positions[id0] += deltaLambda * w0 * gradC;
        Positions[id1] -= deltaLambda * w1 * gradC;
    }
}

/// <summary>
/// 针对解析碰撞体（Sphere / Box）的位置修正 + 可选切向摩擦。
/// 对应旧版：ClothAnalyticalCollisionJob（无摩擦） / SoftBodyAnalyticalCollisionJob（含摩擦）。
/// EnableFriction=false 时等价于旧的 Cloth 版本（不触碰 PrevPositions）。
/// EnableFriction=true  时等价于旧的 SoftBody 版本（通过 PrevPositions 实现切向摩擦）。
/// </summary>
[BurstCompile]
public struct XPBDAnalyticalCollisionJob : IJob
{
    public NativeArray<float3> Positions;
    public NativeArray<float3> PrevPositions;
    [ReadOnly] public NativeArray<float> InvMasses;
    [ReadOnly] public NativeArray<AnalyticalColliderData> Colliders;
    public float ParticleRadius;
    public float Friction;
    public int NumParticles;
    public bool EnableFriction;

    public void Execute()
    {
        for (int i = 0; i < NumParticles; i++)
        {
            if (InvMasses[i] == 0f) continue;

            float3 pos = Positions[i];

            for (int c = 0; c < Colliders.Length; c++)
            {
                var col = Colliders[c];

                float3 newPos;
                float3 pushDir;
                float pushDist;
                bool hit;

                switch (col.Type)
                {
                    case AnalyticalColliderType.Sphere:
                        hit = ResolveSphereCollision(pos, col.Center, col.Radius, ParticleRadius,
                            out newPos, out pushDir, out pushDist);
                        break;
                    case AnalyticalColliderType.Box:
                        hit = ResolveBoxCollision(pos, col, ParticleRadius,
                            out newPos, out pushDir, out pushDist);
                        break;
                    default:
                        hit = false;
                        newPos = pos;
                        pushDir = float3.zero;
                        pushDist = 0f;
                        break;
                }

                if (hit)
                {
                    pos = newPos;

                    if (EnableFriction)
                    {
                        float3 prev = PrevPositions[i];
                        float3 disp = pos - prev;
                        // 去除法向分量（沿 pushDir），只对切向分量做衰减
                        float3 tangent = disp - math.dot(disp, pushDir) * pushDir;
                        float tangentLen = math.length(tangent);

                        if (tangentLen < 1e-4f)
                        {
                            // 切向几乎为零：直接冻结，避免自旋抖动
                            PrevPositions[i] = pos;
                        }
                        else
                        {
                            float k = math.saturate(Friction);
                            PrevPositions[i] += tangent * k;
                        }
                    }
                }
            }

            Positions[i] = pos;
        }
    }

    static bool ResolveSphereCollision(float3 particlePos, float3 sphereCenter, float sphereRadius, float particleRadius,
        out float3 newPos, out float3 pushDir, out float pushDist)
    {
        float3 diff = particlePos - sphereCenter;
        float dist = math.length(diff);
        float minDist = sphereRadius + particleRadius;

        if (dist < minDist && dist > 1e-8f)
        {
            float3 dir = diff / dist;
            newPos = sphereCenter + dir * minDist;
            pushDir = dir;
            pushDist = minDist - dist;
            return true;
        }
        else if (dist <= 1e-8f)
        {
            pushDir = new float3(0, 1, 0);
            newPos = sphereCenter + pushDir * minDist;
            pushDist = minDist;
            return true;
        }

        newPos = particlePos;
        pushDir = float3.zero;
        pushDist = 0f;
        return false;
    }

    static bool ResolveBoxCollision(float3 particlePos, AnalyticalColliderData box, float particleRadius,
        out float3 newPos, out float3 pushDir, out float pushDist)
    {
        float3 localPos = math.mul(box.InvRotation, particlePos - box.Center);
        float3 expandedHalf = box.HalfExtents + particleRadius;

        if (math.abs(localPos.x) < expandedHalf.x &&
            math.abs(localPos.y) < expandedHalf.y &&
            math.abs(localPos.z) < expandedHalf.z)
        {
            float3 penetration = expandedHalf - math.abs(localPos);
            float3 sign = math.sign(localPos);
            float3 localNormal;
            float dist;

            if (penetration.x <= penetration.y && penetration.x <= penetration.z)
            {
                localPos.x = sign.x * expandedHalf.x;
                localNormal = new float3(sign.x, 0, 0);
                dist = penetration.x;
            }
            else if (penetration.y <= penetration.x && penetration.y <= penetration.z)
            {
                localPos.y = sign.y * expandedHalf.y;
                localNormal = new float3(0, sign.y, 0);
                dist = penetration.y;
            }
            else
            {
                localPos.z = sign.z * expandedHalf.z;
                localNormal = new float3(0, 0, sign.z);
                dist = penetration.z;
            }

            newPos = math.mul(box.Rotation, localPos) + box.Center;
            pushDir = math.mul(box.Rotation, localNormal);
            pushDist = dist;
            return true;
        }

        newPos = particlePos;
        pushDir = float3.zero;
        pushDist = 0f;
        return false;
    }
}

/// <summary>
/// 阻尼 + 速度回写。
/// 对应旧版：ClothPostSolveJob / SoftBodyPostSolveJob。
/// </summary>
[BurstCompile]
public struct XPBDPostSolveJob : IJobParallelFor
{
    public NativeArray<float3> Positions;
    [ReadOnly] public NativeArray<float3> PrevPositions;
    public NativeArray<float3> Velocities;
    [ReadOnly] public NativeArray<float> InvMasses;
    public float OneOverDt;
    public float Dt;
    public float Damping;

    public void Execute(int i)
    {
        if (InvMasses[i] == 0f) return;

        if (Damping > 0f)
        {
            float3 displacement = Positions[i] - PrevPositions[i];
            Positions[i] = PrevPositions[i] + displacement * (1f - Damping);
        }

        Velocities[i] = (Positions[i] - PrevPositions[i]) * OneOverDt;
    }
}
