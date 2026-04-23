using Unity.Burst;
using Unity.Collections;
using Unity.Jobs;
using Unity.Mathematics;

// ============================================================
// SoftBody PreSolve Job - 重力积分（可并行）
// ============================================================
[BurstCompile]
public struct SoftBodyPreSolveJob : IJobParallelFor
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
// SoftBody Lambda重置 Job（距离 + 体积）
// ============================================================
[BurstCompile]
public struct SoftBodyResetDistanceLambdaJob : IJobParallelFor
{
    public NativeArray<float> Lambdas;

    public void Execute(int i)
    {
        Lambdas[i] = 0f;
    }
}

[BurstCompile]
public struct SoftBodyResetVolumeLambdaJob : IJobParallelFor
{
    public NativeArray<float> Lambdas;

    public void Execute(int i)
    {
        Lambdas[i] = 0f;
    }
}

// ============================================================
// SoftBody DistanceConstraint Job（顺序依赖，使用IJob + Burst）
// 与Cloth的距离约束完全一致
// 直接消费 SoftBodyEdge 结构体 Buffer（零拷贝），避免每帧拆包成 int/float 数组。
// ============================================================
[BurstCompile]
public struct SoftBodyDistanceConstraintJob : IJob
{
    public NativeArray<float3> Positions;
    [ReadOnly] public NativeArray<float> InvMasses;
    [ReadOnly] public NativeArray<SoftBodyEdge> Edges;
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

// ============================================================
// SoftBody VolumeConstraint Job（四面体体积保持约束）
// XPBD体积约束：C = V - V_rest
// 梯度方向为对面三角形的法向量（面积加权）
// ============================================================
[BurstCompile]
public struct SoftBodyVolumeConstraintJob : IJob
{
    public NativeArray<float3> Positions;
    [ReadOnly] public NativeArray<float> InvMasses;
    [ReadOnly] public NativeArray<Tetrahedron> Tets;
    [ReadOnly] public NativeArray<TetrahedronRestVolume> RestVolumes;
    public NativeArray<float> Lambdas;
    public float Stiffness; // compliance
    public float Dt;

    public void Execute()
    {
        float alpha = Stiffness / (Dt * Dt);

        int tetCount = Tets.Length;
        for (int t = 0; t < tetCount; t++)
        {
            var tet = Tets[t];
            int i0 = tet.I0;
            int i1 = tet.I1;
            int i2 = tet.I2;
            int i3 = tet.I3;

            float3 p0 = Positions[i0];
            float3 p1 = Positions[i1];
            float3 p2 = Positions[i2];
            float3 p3 = Positions[i3];

            // 当前体积 V = dot(p1-p0, cross(p2-p0, p3-p0)) / 6
            float3 d1 = p1 - p0;
            float3 d2 = p2 - p0;
            float3 d3 = p3 - p0;
            float vol = math.dot(d1, math.cross(d2, d3)) / 6f;

            float restVol = RestVolumes[t].Value;
            float C = vol - restVol;

            // 梯度：对每个顶点的位置求导
            // grad_p0 = -cross(p2-p1, p3-p1) / 6
            // grad_p1 =  cross(p2-p0, p3-p0) / 6
            // grad_p2 =  cross(p3-p0, p1-p0) / 6  
            // grad_p3 =  cross(p1-p0, p2-p0) / 6
            float3 grad0 = -math.cross(p2 - p1, p3 - p1) / 6f;
            float3 grad1 =  math.cross(p2 - p0, p3 - p0) / 6f;
            float3 grad2 =  math.cross(p3 - p0, p1 - p0) / 6f;
            float3 grad3 =  math.cross(p1 - p0, p2 - p0) / 6f;

            float w0 = InvMasses[i0];
            float w1 = InvMasses[i1];
            float w2 = InvMasses[i2];
            float w3 = InvMasses[i3];

            float wSum = w0 * math.lengthsq(grad0) +
                         w1 * math.lengthsq(grad1) +
                         w2 * math.lengthsq(grad2) +
                         w3 * math.lengthsq(grad3);

            if (wSum < 1e-12f) continue;

            float deltaLambda = -(C + Lambdas[t] * alpha) / (wSum + alpha);
            Lambdas[t] += deltaLambda;

            Positions[i0] += deltaLambda * w0 * grad0;
            Positions[i1] += deltaLambda * w1 * grad1;
            Positions[i2] += deltaLambda * w2 * grad2;
            Positions[i3] += deltaLambda * w3 * grad3;
        }
    }
}

// ============================================================
// SoftBody 解析碰撞约束 Job（Burst编译，在SubStep内运行）
// 支持球体和Box碰撞体
// 增加了摩擦处理：通过修正 PrevPositions 来耗散切向速度，
// 避免软体球在墙面/地面上无阻力地持续滑动或自旋。
// ============================================================
[BurstCompile]
public struct SoftBodyAnalyticalCollisionJob : IJob
{
    public NativeArray<float3> Positions;
    public NativeArray<float3> PrevPositions;
    [ReadOnly] public NativeArray<float> InvMasses;
    [ReadOnly] public NativeArray<AnalyticalColliderData> Colliders;
    public float ParticleRadius;
    public float Friction; // 0 = 无摩擦，1 = 完全吸附（切向速度全部耗散）
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

                float3 newPos = pos;
                float3 pushDir = float3.zero;
                float pushDist = 0f;
                bool hit = false;

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
                }

                if (hit)
                {
                    pos = newPos;

                    // === 摩擦：把接触点的切向位移按 Friction 比例拉近 prevPos ===
                    // 这等价于让 PostSolve 中 v = (pos - prev)/dt 的切向分量被削减，
                    // 从而耗散沿墙面滑动/旋转的动能。
                    float3 prev = PrevPositions[i];
                    float3 disp = pos - prev;
                    // 去除法向分量（沿 pushDir），只对切向分量做衰减
                    float3 tangent = disp - math.dot(disp, pushDir) * pushDir;
                    float tangentLen = math.length(tangent);

                    // 切向运动很小时直接完全停止（避免抖动、自旋）
                    if (tangentLen < 1e-4f)
                    {
                        PrevPositions[i] = pos;
                    }
                    else
                    {
                        // 把 prev 沿切向朝 pos 拉近 Friction 比例
                        // 系数也会随穿透深度略微放大（穿透越深，摩擦越强）
                        float k = math.saturate(Friction);
                        PrevPositions[i] += tangent * k;
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

// ============================================================
// SoftBody PostSolve Job - 速度更新 + 阻尼（可并行）
// ============================================================
[BurstCompile]
public struct SoftBodyPostSolveJob : IJobParallelFor
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
