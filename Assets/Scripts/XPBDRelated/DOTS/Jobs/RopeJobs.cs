using Unity.Burst;
using Unity.Collections;
using Unity.Jobs;
using Unity.Mathematics;

[BurstCompile]
public struct RopePreSolvePointJob : IJobParallelFor
{
    public NativeArray<float3> PointPos;
    public NativeArray<float3> PrevPos;
    public NativeArray<float3> Vel;
    [ReadOnly] public NativeArray<float> PointInvMass;

    public float3 Gravity;
    public float GravityFactor;
    public float Dt;

    public void Execute(int i)
    {
        if (PointInvMass[i] == 0f)
        {
            PrevPos[i] = PointPos[i];
            return;
        }

        // NaN安全检查：如果位置或速度已经是NaN，重置为安全值
        if (math.any(math.isnan(PointPos[i])) || math.any(math.isnan(Vel[i])))
        {
            Vel[i] = float3.zero;
            PrevPos[i] = PointPos[i];
            return;
        }

        float3 vel = Vel[i];
        vel += Gravity * GravityFactor * Dt;

        PrevPos[i] = PointPos[i];

        // ClampMagnitude
        float maxSpeed = math.sqrt(2.5f * 2.5f * 3f);
        float speed = math.length(vel);
        if (speed > maxSpeed)
            vel = math.normalize(vel) * maxSpeed;

        Vel[i] = vel;
        PointPos[i] += vel * Dt;
    }
}

[BurstCompile]
public struct RopePreSolveGhostJob : IJobParallelFor
{
    public NativeArray<float3> GhostPos;
    public NativeArray<float3> GhostPrev;
    [ReadOnly] public NativeArray<float3> GhostVels;
    public float Dt;

    public void Execute(int i)
    {
        GhostPrev[i] = GhostPos[i];
        GhostPos[i] += GhostVels[i] * Dt;
    }
}

[BurstCompile]
public struct RopeResetLambdaJob : IJobParallelFor
{
    public NativeArray<float> EdgeLambda0;
    public NativeArray<float> EdgeLambda1;
    public NativeArray<float> EdgeLambda2;

    public void Execute(int i)
    {
        EdgeLambda0[i] = 0f;
        EdgeLambda1[i] = 0f;
        EdgeLambda2[i] = 0f;
    }
}

[BurstCompile]
public struct RopeResetBendLambdaJob : IJobParallelFor
{
    public NativeArray<float3> BendLambdas;

    public void Execute(int i)
    {
        BendLambdas[i] = float3.zero;
    }
}

[BurstCompile]
public struct RopeEdgeConstraintJob : IJob
{
    public NativeArray<float3> PointPos;
    public NativeArray<float3> GhostPos;
    [ReadOnly] public NativeArray<float> PointInvMass;
    [ReadOnly] public NativeArray<float> GhostInvMass;
    [ReadOnly] public NativeArray<float> RestLengths;
    public NativeArray<float> Lambda0;
    public NativeArray<float> Lambda1;
    public NativeArray<float> Lambda2;
    public float Stiffness;
    public float GhostDistance;
    public float Dt;
    public int NumPoints;
    [ReadOnly] public NativeArray<byte> SegmentContactFlags;

    public void Execute()
    {
        const float EPSILON = 1e-6f;
        float alpha = Stiffness / (Dt * Dt);

        int numEdges = NumPoints - 1;
        for (int step = 0; step < numEdges; step++)
        {
            int i;
            int half = step >> 1;                 // step / 2
            if ((step & 1) == 0) i = half;         // 偶数步：从左端向中间
            else i = numEdges - 1 - half;          // 奇数步：从右端向中间

            // NaN安全检查：阻断NaN传播
            if (math.any(math.isnan(PointPos[i])) || math.any(math.isnan(PointPos[i + 1])) ||
                math.any(math.isnan(GhostPos[i])))
                continue;

            // === 约束1：相邻点距离约束 ===
            float3 dir = PointPos[i] - PointPos[i + 1];
            float len = math.length(dir);
            float wSum = PointInvMass[i] + PointInvMass[i + 1];

            if (len > EPSILON && wSum > EPSILON)
            {
                float deltaLambda = -((len - RestLengths[i]) + alpha * Lambda0[i]) / (wSum + alpha);
                float3 dP = deltaLambda * (dir / len);
                PointPos[i] += dP * PointInvMass[i];
                PointPos[i + 1] -= dP * PointInvMass[i + 1];
                Lambda0[i] += deltaLambda;
            }

            // 段一旦处于接触中，只保留长度约束；
            // 跳过 ghost 相关约束，避免辅助框架持续把中心线拉回碰撞体表面。
            bool segmentInContact = i < SegmentContactFlags.Length && SegmentContactFlags[i] != 0;
            if (segmentInContact)
                continue;

            // === 约束2：Ghost点共面约束 ===
            float3 pm = 0.5f * (PointPos[i] + PointPos[i + 1]);
            float3 p0p2 = PointPos[i] - GhostPos[i];
            float3 p2p1 = GhostPos[i] - PointPos[i + 1];
            float3 p1p0 = PointPos[i + 1] - PointPos[i];
            float3 p2pm = GhostPos[i] - pm;

            wSum = PointInvMass[i] * math.lengthsq(p0p2) +
                   PointInvMass[i + 1] * math.lengthsq(p2p1) +
                   GhostInvMass[i] * math.lengthsq(p1p0);

            if (wSum > EPSILON)
            {
                float p0p2_len = math.length(p0p2);
                float p2p1_len = math.length(p2p1);
                float p1p0_len = math.length(p1p0);

                // 防止normalize零向量产生NaN
                if (p0p2_len > EPSILON && p2p1_len > EPSILON && p1p0_len > EPSILON)
                {
                    float deltaLambda = -(math.dot(p2pm, p1p0) + alpha * Lambda1[i]) / (wSum + alpha);
                    PointPos[i] += (p0p2 / p0p2_len) * deltaLambda * PointInvMass[i];
                    PointPos[i + 1] += (p2p1 / p2p1_len) * deltaLambda * PointInvMass[i + 1];
                    GhostPos[i] += (p1p0 / p1p0_len) * deltaLambda * GhostInvMass[i];
                    Lambda1[i] += deltaLambda;
                }
            }

            // === 约束3：Ghost点距离约束 ===
            wSum = 0.25f * PointInvMass[i] + 0.25f * PointInvMass[i + 1] + 1.0f * GhostInvMass[i];

            if (wSum > EPSILON)
            {
                pm = 0.5f * (PointPos[i] + PointPos[i + 1]);
                p2pm = GhostPos[i] - pm;

                float p2pm_mag = math.length(p2pm);

                // 防止除零产生NaN
                if (p2pm_mag > EPSILON)
                {
                    float3 p2pm_dir = p2pm * (1.0f / p2pm_mag);

                    float deltaLambda = -(p2pm_mag - GhostDistance + alpha * Lambda2[i]) / (wSum + alpha);

                    PointPos[i] -= 0.5f * PointInvMass[i] * deltaLambda * p2pm_dir;
                    PointPos[i + 1] -= 0.5f * PointInvMass[i + 1] * deltaLambda * p2pm_dir;
                    GhostPos[i] += 1.0f * GhostInvMass[i] * deltaLambda * p2pm_dir;
                    Lambda2[i] += deltaLambda;
                }
            }
        }
    }
}

[BurstCompile]
public struct RopeBendTwistConstraintJob : IJob
{
    // 排列组合常量
    // permutation[0] = {0,2,1}, permutation[1] = {1,0,2}, permutation[2] = {2,1,0}

    public NativeArray<float3> PointPos;
    public NativeArray<float3> GhostPos;
    [ReadOnly] public NativeArray<float> PointInvMass;
    [ReadOnly] public NativeArray<float> GhostInvMass;
    [ReadOnly] public NativeArray<float3> InitDarboux;
    [ReadOnly] public NativeArray<float> RestLengths;
    public NativeArray<float3> Lambdas;
    public float BendTwistKs;
    public float Stiffness;
    public float Dt;
    public int NumPoints;
    [ReadOnly] public NativeArray<byte> SegmentContactFlags;

    public void Execute()
    {
        float alpha = Stiffness / (Dt * Dt);

        int numBend = NumPoints - 2;
        for (int step = 0; step < numBend; step++)
        {
            int idx;
            int half = step >> 1;
            if ((step & 1) == 0) idx = half;
            else idx = numBend - 1 - half;

            bool leftSegmentInContact = idx < SegmentContactFlags.Length && SegmentContactFlags[idx] != 0;
            bool rightSegmentInContact = (idx + 1) < SegmentContactFlags.Length && SegmentContactFlags[idx + 1] != 0;
            if (leftSegmentInContact || rightSegmentInContact)
                continue;

            float halfLength = RestLengths[0] / 2f;
            if (halfLength < 1e-8f) continue;

            // 检查输入是否包含NaN（防止NaN传播）
            if (math.any(math.isnan(PointPos[idx])) || math.any(math.isnan(PointPos[idx + 1])) ||
                math.any(math.isnan(PointPos[idx + 2])) || math.any(math.isnan(GhostPos[idx])) ||
                math.any(math.isnan(GhostPos[idx + 1])))
                continue;

            // 检查退化情况：相邻点重合或Ghost点与线段共线
            float len01 = math.length(PointPos[idx + 1] - PointPos[idx]);
            float len12 = math.length(PointPos[idx + 2] - PointPos[idx + 1]);
            if (len01 < 1e-6f || len12 < 1e-6f) continue;

            // 计算材料坐标系
            float3x3 a = ComputeMaterialFrame(PointPos[idx], PointPos[idx + 1], GhostPos[idx]);
            float3x3 b = ComputeMaterialFrame(PointPos[idx + 1], PointPos[idx + 2], GhostPos[idx + 1]);

            if (math.lengthsq(a.c0) < 0.5f || math.lengthsq(b.c0) < 0.5f) continue;

            // 计算Darboux向量
            float3 darboux = ComputeDarbouxVector(a, b, halfLength);
            if (math.any(math.isnan(darboux))) continue;

            float3x3 da0p0, da0p1, da0p2;
            float3x3 da1p0, da1p1, da1p2;
            float3x3 da2p0, da2p1, da2p2;
            ComputeMaterialFrameDerivative(PointPos[idx], PointPos[idx + 1], GhostPos[idx], a,
                out da0p0, out da0p1, out da0p2,
                out da1p0, out da1p1, out da1p2,
                out da2p0, out da2p1, out da2p2);

            float3x3 db0p0, db0p1, db0p2;
            float3x3 db1p0, db1p1, db1p2;
            float3x3 db2p0, db2p1, db2p2;
            ComputeMaterialFrameDerivative(PointPos[idx + 1], PointPos[idx + 2], GhostPos[idx + 1], b,
                out db0p0, out db0p1, out db0p2,
                out db1p0, out db1p1, out db1p2,
                out db2p0, out db2p1, out db2p2);

            // 计算Darboux梯度
            float3x3 omega_pa, omega_pb, omega_pc, omega_pd, omega_pe;
            ComputeDarbouxGradient(darboux, halfLength, a, b,
                da0p0, da0p1, da0p2, da1p0, da1p1, da1p2, da2p0, da2p1, da2p2,
                db0p0, db0p1, db0p2, db1p0, db1p1, db1p2, db2p0, db2p1, db2p2,
                BendTwistKs,
                out omega_pa, out omega_pb, out omega_pc, out omega_pd, out omega_pe);

            // 计算约束值 C
            float3 C = new float3(
                (darboux.x - InitDarboux[idx].x),
                (darboux.y - InitDarboux[idx].y),
                (darboux.z - InitDarboux[idx].z)
            ) * BendTwistKs;

            // 计算factor_matrix
            float3x3 factor_matrix = float3x3.zero;

            float5 invMasses = new float5(
                PointInvMass[idx], PointInvMass[idx + 1], PointInvMass[idx + 2],
                GhostInvMass[idx], GhostInvMass[idx + 1]
            );

            // 对5个jacobian矩阵累加
            AccumulateFactorMatrix(ref factor_matrix, omega_pa, invMasses.x);
            AccumulateFactorMatrix(ref factor_matrix, omega_pb, invMasses.y);
            AccumulateFactorMatrix(ref factor_matrix, omega_pc, invMasses.z);
            AccumulateFactorMatrix(ref factor_matrix, omega_pd, invMasses.w);
            AccumulateFactorMatrix(ref factor_matrix, omega_pe, invMasses.v);

            // 加上alpha项
            factor_matrix.c0 *= alpha;
            factor_matrix.c1 *= alpha;
            factor_matrix.c2 *= alpha;

            // 求逆
            float3x3 inv_factor = Inverse3x3(factor_matrix);

            // deltaLambda = inv_factor * -(C + alpha * lambdas[idx])
            float3 rhs = -(C + alpha * Lambdas[idx]);
            float3 deltaLambda = math.mul(inv_factor, rhs);

            // NaN安全检查：如果deltaLambda包含NaN或Inf，跳过此约束
            if (math.any(math.isnan(deltaLambda)) || math.any(math.isinf(deltaLambda)))
                continue;

            Lambdas[idx] += deltaLambda;

            // 计算位置修正
            float3 dp0 = math.mul(ScaleColumns(omega_pa, invMasses.x), deltaLambda);
            float3 dp1 = math.mul(ScaleColumns(omega_pb, invMasses.y), deltaLambda);
            float3 dp2 = math.mul(ScaleColumns(omega_pc, invMasses.z), deltaLambda);
            float3 dp3 = math.mul(ScaleColumns(omega_pd, invMasses.w), deltaLambda);
            float3 dp4 = math.mul(ScaleColumns(omega_pe, invMasses.v), deltaLambda);

            PointPos[idx] += dp0;
            PointPos[idx + 1] += dp1;
            PointPos[idx + 2] += dp2;
            GhostPos[idx] += dp3;
            GhostPos[idx + 1] += dp4;
        }
    }

    // 辅助结构体，存储5个float
    struct float5
    {
        public float x, y, z, w, v;
        public float5(float x, float y, float z, float w, float v)
        {
            this.x = x; this.y = y; this.z = z; this.w = w; this.v = v;
        }
    }

    static float3x3 ScaleColumns(float3x3 m, float s)
    {
        return new float3x3(m.c0 * s, m.c1 * s, m.c2 * s);
    }

    static void AccumulateFactorMatrix(ref float3x3 factor, float3x3 jacobian, float invMass)
    {
        float3x3 jt = math.transpose(jacobian);
        float3x3 jtj = math.mul(jt, jacobian);
        float3x3 scaled = ScaleColumns(jtj, invMass);
        factor = new float3x3(
            factor.c0 + scaled.c0,
            factor.c1 + scaled.c1,
            factor.c2 + scaled.c2
        );
    }

    static float3x3 ComputeMaterialFrame(float3 p1, float3 p2, float3 pg)
    {
        float3 diff = p2 - p1;
        float diffLen = math.length(diff);
        if (diffLen < 1e-8f) return float3x3.identity;

        float3 d3 = diff / diffLen;
        float3 crossVal = math.cross(d3, pg - p1);
        float crossLen = math.length(crossVal);
        if (crossLen < 1e-8f) return float3x3.identity;

        float3 d2 = crossVal / crossLen;
        float3 d1 = math.cross(d2, d3);
        // 列存储：c0=d1, c1=d2, c2=d3
        return new float3x3(d1, d2, d3);
    }

    static float3 ComputeDarbouxVector(float3x3 dA, float3x3 dB, float halfLength)
    {
        float factor = 1.0f + math.dot(dA.c0, dB.c0) + math.dot(dA.c1, dB.c1) + math.dot(dA.c2, dB.c2);

        // 防止除零（当两个frame几乎反向时factor接近0）
        if (math.abs(factor) < 1e-8f) return float3.zero;

        factor = 2.0f / (halfLength * factor);

        // permutation: {0,2,1}, {1,0,2}, {2,1,0}
        float3 darboux;
        darboux.x = math.dot(dA.c2, dB.c1) - math.dot(dA.c1, dB.c2); // i=0, j=2, k=1
        darboux.y = math.dot(dA.c0, dB.c2) - math.dot(dA.c2, dB.c0); // i=1, j=0, k=2
        darboux.z = math.dot(dA.c1, dB.c0) - math.dot(dA.c0, dB.c1); // i=2, j=1, k=0

        return darboux * factor;
    }

    static void ComputeMaterialFrameDerivative(float3 p0, float3 p1, float3 p2, float3x3 d,
        out float3x3 d1p0, out float3x3 d1p1, out float3x3 d1p2,
        out float3x3 d2p0, out float3x3 d2p1, out float3x3 d2p2,
        out float3x3 d3p0, out float3x3 d3p1, out float3x3 d3p2)
    {
        // d3pi
        float3 p01 = p1 - p0;
        float length_p01 = math.length(p01);

        // 防止除零
        if (length_p01 < 1e-8f) length_p01 = 1e-8f;

        // d3p0 = (d3 * d3^T - I) / length_p01
        float3x3 d3outer = OuterProduct(d.c2, d.c2);
        d3p0 = new float3x3(
            (d3outer.c0 - new float3(1, 0, 0)) / length_p01,
            (d3outer.c1 - new float3(0, 1, 0)) / length_p01,
            (d3outer.c2 - new float3(0, 0, 1)) / length_p01
        );
        d3p1 = new float3x3(-d3p0.c0, -d3p0.c1, -d3p0.c2);
        d3p2 = float3x3.zero;

        // d2pi
        float3 p02 = p2 - p0;
        float3 p01_cross_p02 = math.cross(p01, p02);
        float length_cross = math.length(p01_cross_p02);

        // 防止除零
        if (length_cross < 1e-8f) length_cross = 1e-8f;

        float3x3 d2outer = OuterProduct(d.c1, d.c1);
        float3x3 mat = new float3x3(
            (d2outer.c0 - new float3(1, 0, 0)) * (-1.0f / length_cross),
            (d2outer.c1 - new float3(0, 1, 0)) * (-1.0f / length_cross),
            (d2outer.c2 - new float3(0, 0, 1)) * (-1.0f / length_cross)
        );

        d2p0 = math.mul(mat, CrossMatrix(p2 - p1));
        d2p1 = math.mul(mat, CrossMatrix(p0 - p2));
        d2p2 = math.mul(mat, CrossMatrix(p1 - p0));

        // d1pi = cross(d2, d3) 的导数
        float3x3 crossD3 = CrossMatrix(d.c2);
        float3x3 crossD2 = CrossMatrix(d.c1);

        d1p0 = math.mul(crossD2, d3p0) - math.mul(crossD3, d2p0);
        d1p1 = math.mul(crossD2, d3p1) - math.mul(crossD3, d2p1);
        // d3p2 = 0, 所以 crossD2 * d3p2 = 0
        d1p2 = new float3x3(-math.mul(crossD3, d2p2).c0, -math.mul(crossD3, d2p2).c1, -math.mul(crossD3, d2p2).c2);
    }

    static void ComputeDarbouxGradient(
        float3 darboux_vector, float halfLength,
        float3x3 da, float3x3 db,
        float3x3 da0p0, float3x3 da0p1, float3x3 da0p2,
        float3x3 da1p0, float3x3 da1p1, float3x3 da1p2,
        float3x3 da2p0, float3x3 da2p1, float3x3 da2p2,
        float3x3 db0p0, float3x3 db0p1, float3x3 db0p2,
        float3x3 db1p0, float3x3 db1p1, float3x3 db1p2,
        float3x3 db2p0, float3x3 db2p1, float3x3 db2p2,
        float bendKs,
        out float3x3 omega_pa, out float3x3 omega_pb, out float3x3 omega_pc,
        out float3x3 omega_pd, out float3x3 omega_pe)
    {
        omega_pa = float3x3.zero;
        omega_pb = float3x3.zero;
        omega_pc = float3x3.zero;
        omega_pd = float3x3.zero;
        omega_pe = float3x3.zero;

        float x = 1.0f + math.dot(da.c0, db.c0) + math.dot(da.c1, db.c1) + math.dot(da.c2, db.c2);
        float denom = halfLength * x;
        // 防止除零
        if (math.abs(denom) < 1e-8f) return;
        x = 2.0f / denom;

        // 我们需要按 j,k 索引访问

        for (int c = 0; c < 3; c++)
        {
            int ii, jj, kk;
            GetPermutation(c, out ii, out jj, out kk);

            // === omega_pa (point a = p0, 只有da的导数对p0) ===
            {
                float3x3 dajp0 = GetFrameDerivative(jj, 0, da0p0, da1p0, da2p0);
                float3x3 dakp0 = GetFrameDerivative(kk, 0, da0p0, da1p0, da2p0);

                float3 term1 = math.mul(math.transpose(dajp0), GetColumn(db, kk))
                             - math.mul(math.transpose(dakp0), GetColumn(db, jj));

                float3 term2 = float3.zero;
                for (int n = 0; n < 3; n++)
                {
                    float3x3 danp0 = GetFrameDerivative(n, 0, da0p0, da1p0, da2p0);
                    term2 += math.mul(math.transpose(danp0), GetColumn(db, n));
                }

                float3 col = (term1 - 0.5f * darboux_vector[ii] * halfLength * term2) * x * bendKs;
                SetColumn(ref omega_pa, ii, col);
            }

            // === omega_pb (point b = p1, 有da对p1和db对p0的导数) ===
            {
                float3x3 dajp1 = GetFrameDerivative(jj, 1, da0p1, da1p1, da2p1);
                float3x3 dakp1 = GetFrameDerivative(kk, 1, da0p1, da1p1, da2p1);
                float3x3 dbjp0 = GetFrameDerivative(jj, 0, db0p0, db1p0, db2p0);
                float3x3 dbkp0 = GetFrameDerivative(kk, 0, db0p0, db1p0, db2p0);

                float3 term1 = math.mul(math.transpose(dajp1), GetColumn(db, kk))
                             - math.mul(math.transpose(dakp1), GetColumn(db, jj))
                             - math.mul(math.transpose(dbjp0), GetColumn(da, kk))
                             + math.mul(math.transpose(dbkp0), GetColumn(da, jj));

                float3 term2 = float3.zero;
                for (int n = 0; n < 3; n++)
                {
                    float3x3 danp1 = GetFrameDerivative(n, 1, da0p1, da1p1, da2p1);
                    float3x3 dbnp0 = GetFrameDerivative(n, 0, db0p0, db1p0, db2p0);
                    term2 += math.mul(math.transpose(danp1), GetColumn(db, n));
                    term2 += math.mul(math.transpose(dbnp0), GetColumn(da, n));
                }

                float3 col = (term1 - 0.5f * darboux_vector[ii] * halfLength * term2) * x * bendKs;
                SetColumn(ref omega_pb, ii, col);
            }

            // === omega_pc (point c = p2, 只有db对p1的导数) ===
            {
                float3x3 dbjp1 = GetFrameDerivative(jj, 1, db0p1, db1p1, db2p1);
                float3x3 dbkp1 = GetFrameDerivative(kk, 1, db0p1, db1p1, db2p1);

                float3 term1 = math.mul(math.transpose(dbjp1), GetColumn(da, kk))
                             - math.mul(math.transpose(dbkp1), GetColumn(da, jj));

                float3 term2 = float3.zero;
                for (int n = 0; n < 3; n++)
                {
                    float3x3 dbnp1 = GetFrameDerivative(n, 1, db0p1, db1p1, db2p1);
                    term2 += math.mul(math.transpose(dbnp1), GetColumn(da, n));
                }

                float3 col = (term1 + 0.5f * darboux_vector[ii] * halfLength * term2) * (-x) * bendKs;
                SetColumn(ref omega_pc, ii, col);
            }

            // === omega_pd (ghost point d = ghost[idx], 只有da对p2的导数) ===
            {
                float3x3 dajp2 = GetFrameDerivative(jj, 2, da0p2, da1p2, da2p2);
                float3x3 dakp2 = GetFrameDerivative(kk, 2, da0p2, da1p2, da2p2);

                float3 term1 = math.mul(math.transpose(dajp2), GetColumn(db, kk))
                             - math.mul(math.transpose(dakp2), GetColumn(db, jj));

                float3 term2 = float3.zero;
                for (int n = 0; n < 3; n++)
                {
                    float3x3 danp2 = GetFrameDerivative(n, 2, da0p2, da1p2, da2p2);
                    term2 += math.mul(math.transpose(danp2), GetColumn(db, n));
                }

                float3 col = (term1 - 0.5f * darboux_vector[ii] * halfLength * term2) * x * bendKs;
                SetColumn(ref omega_pd, ii, col);
            }

            // === omega_pe (ghost point e = ghost[idx+1], 只有db对p2的导数) ===
            {
                float3x3 dbjp2 = GetFrameDerivative(jj, 2, db0p2, db1p2, db2p2);
                float3x3 dbkp2 = GetFrameDerivative(kk, 2, db0p2, db1p2, db2p2);

                float3 term1 = math.mul(math.transpose(dbjp2), GetColumn(da, kk))
                             - math.mul(math.transpose(dbkp2), GetColumn(da, jj));

                float3 term2 = float3.zero;
                for (int n = 0; n < 3; n++)
                {
                    float3x3 dbnp2 = GetFrameDerivative(n, 2, db0p2, db1p2, db2p2);
                    term2 += math.mul(math.transpose(dbnp2), GetColumn(da, n));
                }

                float3 col = (term1 + 0.5f * darboux_vector[ii] * halfLength * term2) * (-x) * bendKs;
                SetColumn(ref omega_pe, ii, col);
            }
        }
    }

    static void GetPermutation(int c, out int i, out int j, out int k)
    {
        // permutation[0] = {0,2,1}, permutation[1] = {1,0,2}, permutation[2] = {2,1,0}
        switch (c)
        {
            case 0: i = 0; j = 2; k = 1; break;
            case 1: i = 1; j = 0; k = 2; break;
            default: i = 2; j = 1; k = 0; break;
        }
    }

    static float3x3 GetFrameDerivative(int axis, int point,
        float3x3 d0px, float3x3 d1px, float3x3 d2px)
    {
        switch (axis)
        {
            case 0: return d0px;
            case 1: return d1px;
            default: return d2px;
        }
    }

    static float3 GetColumn(float3x3 m, int col)
    {
        switch (col)
        {
            case 0: return m.c0;
            case 1: return m.c1;
            default: return m.c2;
        }
    }

    static void SetColumn(ref float3x3 m, int col, float3 value)
    {
        switch (col)
        {
            case 0: m.c0 = value; break;
            case 1: m.c1 = value; break;
            default: m.c2 = value; break;
        }
    }

    static float3x3 OuterProduct(float3 a, float3 b)
    {
        return new float3x3(a * b.x, a * b.y, a * b.z);
    }

    static float3x3 CrossMatrix(float3 v)
    {
        return new float3x3(
            new float3(0, v.z, -v.y),
            new float3(-v.z, 0, v.x),
            new float3(v.y, -v.x, 0)
        );
    }

    static float3x3 Inverse3x3(float3x3 m)
    {
        // 先检查行列式，防止奇异矩阵导致NaN
        float det = math.determinant(new float4x4(
            new float4(m.c0, 0),
            new float4(m.c1, 0),
            new float4(m.c2, 0),
            new float4(0, 0, 0, 1)
        ));
        if (math.abs(det) < 1e-10f) return float3x3.identity;

        // 使用math.inverse对float3x3
        // float3x3没有直接的inverse，需要转为float4x4
        float4x4 m4 = new float4x4(
            new float4(m.c0, 0),
            new float4(m.c1, 0),
            new float4(m.c2, 0),
            new float4(0, 0, 0, 1)
        );
        float4x4 inv4 = math.inverse(m4);
        return new float3x3(inv4.c0.xyz, inv4.c1.xyz, inv4.c2.xyz);
    }
}

[BurstCompile]
public struct RopeAnalyticalCollisionJob : IJob
{
    public NativeArray<float3> PointPos;
    public NativeArray<float3> PrevPos;
    [ReadOnly] public NativeArray<float> PointInvMass;
    [ReadOnly] public NativeArray<AnalyticalColliderData> Colliders;
    public float ParticleRadius;
    public float Friction;
    public int NumPoints;
    public NativeArray<byte> SegmentContactFlags;

    public void Execute()
    {
        float staticThresh = math.max(ParticleRadius * 0.25f, 1e-3f);
        const float MoveEpsilonSq = 1e-12f;

        for (int i = 0; i < SegmentContactFlags.Length; i++)
        {
            SegmentContactFlags[i] = 0;
        }

        for (int i = 0; i < NumPoints - 1; i++)
        {
            float w0 = PointInvMass[i];
            float w1 = PointInvMass[i + 1];
            if (w0 == 0f && w1 == 0f) continue;

            float3 a = PointPos[i];
            float3 b = PointPos[i + 1];
            if (math.any(math.isnan(a)) || math.any(math.isnan(b))) continue;

            float3 prevA = PrevPos[i];
            float3 prevB = PrevPos[i + 1];

            for (int c = 0; c < Colliders.Length; c++)
            {
                var col = Colliders[c];
                float3 correction;
                float3 hitNormal;
                float hitT;
                bool hit = false;

                switch (col.Type)
                {
                    case AnalyticalColliderType.Sphere:
                        hit = ResolveSphereSegmentCollision(a, b, col.Center, col.Radius, ParticleRadius,
                            out correction, out hitNormal, out hitT);
                        break;
                    case AnalyticalColliderType.Box:
                        hit = ResolveBoxSegmentCollision(a, b, col, ParticleRadius,
                            out correction, out hitNormal, out hitT);
                        break;
                    default:
                        correction = float3.zero;
                        hitNormal = float3.zero;
                        hitT = 0.5f;
                        break;
                }

                if (!hit) continue;

                SegmentContactFlags[i] = 1;

                ApplySegmentCorrection(ref a, ref b, w0, w1, hitT, correction,
                    out float3 moveA, out float3 moveB);

                if (Friction > 0f)
                {
                    if (math.lengthsq(moveA) > MoveEpsilonSq)
                        ApplyEndpointFriction(ref prevA, a, hitNormal, staticThresh, Friction);

                    if (math.lengthsq(moveB) > MoveEpsilonSq)
                        ApplyEndpointFriction(ref prevB, b, hitNormal, staticThresh, Friction);
                }
            }

            PointPos[i] = a;
            PointPos[i + 1] = b;
            PrevPos[i] = prevA;
            PrevPos[i + 1] = prevB;
        }
    }

    static void ApplySegmentCorrection(ref float3 a, ref float3 b,
        float w0, float w1, float t, float3 correction,
        out float3 moveA, out float3 moveB)
    {
        if (w0 > 0f && w1 > 0f)
        {
            moveA = correction;
            moveB = correction;
        }
        else if (w0 > 0f)
        {
            moveA = correction;
            moveB = float3.zero;
        }
        else if (w1 > 0f)
        {
            moveA = float3.zero;
            moveB = correction;
        }
        else
        {
            moveA = float3.zero;
            moveB = float3.zero;
        }

        a += moveA;
        b += moveB;
    }

    static void ApplyEndpointFriction(ref float3 prev, float3 pos, float3 normal,
        float staticThresh, float friction)
    {
        float normalLenSq = math.lengthsq(normal);
        if (normalLenSq < 1e-8f) return;

        float3 n = normal * math.rsqrt(normalLenSq);

        float3 disp = pos - prev;
        float normalProj = math.dot(disp, n);
        float3 tangent = disp - normalProj * n;
        float tLen = math.length(tangent);

        if (tLen < staticThresh)
        {
            prev = pos;
        }
        else
        {
            float3 newTangent = tangent * (1f - friction);
            float3 normalComp = normalProj * n;
            prev = pos - (newTangent + normalComp);
        }
    }

    static bool ResolveSphereSegmentCollision(float3 a, float3 b, float3 sphereCenter,
        float sphereRadius, float particleRadius,
        out float3 correction, out float3 normal, out float t)
    {
        const float EPSILON = 1e-8f;

        float3 ab = b - a;
        float abLenSq = math.lengthsq(ab);
        if (abLenSq > EPSILON)
        {
            t = math.clamp(math.dot(sphereCenter - a, ab) / abLenSq, 0f, 1f);
        }
        else
        {
            t = 0f;
        }

        float3 q = a + ab * t;
        float3 diff = q - sphereCenter;
        float dist = math.length(diff);
        float minDist = sphereRadius + particleRadius;

        if (dist < minDist)
        {
            if (dist > EPSILON)
            {
                normal = diff / dist;
            }
            else
            {
                normal = new float3(0f, 1f, 0f);
            }

            correction = normal * (minDist - dist);
            return true;
        }

        correction = float3.zero;
        normal = float3.zero;
        return false;
    }

    static bool ResolveBoxSegmentCollision(float3 aWorld, float3 bWorld, AnalyticalColliderData box,
        float particleRadius, out float3 correction, out float3 normal, out float t)
    {
        const float EPSILON = 1e-8f;

        float3 a = math.mul(box.InvRotation, aWorld - box.Center);
        float3 b = math.mul(box.InvRotation, bWorld - box.Center);
        float3 d = b - a;
        float3 expandedHalf = box.HalfExtents + particleRadius;

        float tEnter = 0f;
        float tExit = 1f;
        int enterAxis = -1;
        float enterSign = 0f;

        for (int axis = 0; axis < 3; axis++)
        {
            float aComp = GetAxis(a, axis);
            float dComp = GetAxis(d, axis);
            float half = GetAxis(expandedHalf, axis);

            if (math.abs(dComp) < EPSILON)
            {
                if (aComp < -half || aComp > half)
                {
                    correction = float3.zero;
                    normal = float3.zero;
                    t = 0f;
                    return false;
                }

                continue;
            }

            float tNear = (-half - aComp) / dComp;
            float tFar = (half - aComp) / dComp;
            float nearSign = -1f;

            if (tNear > tFar)
            {
                float tmp = tNear;
                tNear = tFar;
                tFar = tmp;
                nearSign = 1f;
            }

            if (tNear > tEnter)
            {
                tEnter = tNear;
                enterAxis = axis;
                enterSign = nearSign;
            }

            tExit = math.min(tExit, tFar);

            if (tEnter > tExit)
            {
                correction = float3.zero;
                normal = float3.zero;
                t = 0f;
                return false;
            }
        }

        if (tExit < 0f || tEnter > 1f)
        {
            correction = float3.zero;
            normal = float3.zero;
            t = 0f;
            return false;
        }

        float clampedEnter = math.clamp(tEnter, 0f, 1f);
        float clampedExit = math.clamp(tExit, 0f, 1f);
        if (clampedEnter > clampedExit)
        {
            correction = float3.zero;
            normal = float3.zero;
            t = 0f;
            return false;
        }

        t = 0.5f * (clampedEnter + clampedExit);
        float3 q = a + d * t;

        float3 penetration = expandedHalf - math.abs(q);
        if (penetration.x < 0f || penetration.y < 0f || penetration.z < 0f)
        {
            correction = float3.zero;
            normal = float3.zero;
            t = 0f;
            return false;
        }

        int axisIndex = 0;
        float depth = penetration.x;
        if (penetration.y < depth)
        {
            axisIndex = 1;
            depth = penetration.y;
        }
        if (penetration.z < depth)
        {
            axisIndex = 2;
            depth = penetration.z;
        }

        float sign = math.sign(GetAxis(q, axisIndex));
        if (sign == 0f)
        {
            if (enterAxis == axisIndex && enterSign != 0f)
                sign = enterSign;
            else
            {
                float dirComp = GetAxis(d, axisIndex);
                sign = dirComp >= 0f ? 1f : -1f;
            }
        }

        float3 localNormal = AxisVector(axisIndex) * sign;
        float3 localCorrection = localNormal * depth;

        correction = math.mul(box.Rotation, localCorrection);
        normal = math.mul(box.Rotation, localNormal);
        return true;
    }

    static float GetAxis(float3 v, int axis)
    {
        switch (axis)
        {
            case 0: return v.x;
            case 1: return v.y;
            default: return v.z;
        }
    }

    static float3 AxisVector(int axis)
    {
        switch (axis)
        {
            case 0: return new float3(1f, 0f, 0f);
            case 1: return new float3(0f, 1f, 0f);
            default: return new float3(0f, 0f, 1f);
        }
    }
}

[BurstCompile]
public struct RopeGhostAnalyticalCollisionJob : IJobParallelFor
{
    public NativeArray<float3> GhostPos;
    [ReadOnly] public NativeArray<float> GhostInvMass;
    [ReadOnly] public NativeArray<AnalyticalColliderData> Colliders;
    public float GhostCollisionRadius;

    public void Execute(int i)
    {
        if (GhostInvMass[i] == 0f) return;

        float3 pos = GhostPos[i];
        if (math.any(math.isnan(pos))) return;

        for (int c = 0; c < Colliders.Length; c++)
        {
            var col = Colliders[c];
            switch (col.Type)
            {
                case AnalyticalColliderType.Sphere:
                    pos = ResolveSphereCollision(pos, col.Center, col.Radius, GhostCollisionRadius);
                    break;
                case AnalyticalColliderType.Box:
                    pos = ResolveBoxCollision(pos, col, GhostCollisionRadius);
                    break;
            }
        }

        GhostPos[i] = pos;
    }

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
            particlePos = sphereCenter + new float3(0f, 1f, 0f) * minDist;
        }

        return particlePos;
    }

    static float3 ResolveBoxCollision(float3 particlePos, AnalyticalColliderData box, float particleRadius)
    {
        float3 localPos = math.mul(box.InvRotation, particlePos - box.Center);
        float3 expandedHalf = box.HalfExtents + particleRadius;

        if (math.abs(localPos.x) < expandedHalf.x &&
            math.abs(localPos.y) < expandedHalf.y &&
            math.abs(localPos.z) < expandedHalf.z)
        {
            float3 penetration = expandedHalf - math.abs(localPos);
            int axisIndex = 0;
            float depth = penetration.x;
            if (penetration.y < depth)
            {
                axisIndex = 1;
                depth = penetration.y;
            }
            if (penetration.z < depth)
            {
                axisIndex = 2;
            }

            switch (axisIndex)
            {
                case 0:
                    localPos.x = (localPos.x >= 0f ? 1f : -1f) * expandedHalf.x;
                    break;
                case 1:
                    localPos.y = (localPos.y >= 0f ? 1f : -1f) * expandedHalf.y;
                    break;
                default:
                    localPos.z = (localPos.z >= 0f ? 1f : -1f) * expandedHalf.z;
                    break;
            }

            particlePos = math.mul(box.Rotation, localPos) + box.Center;
        }

        return particlePos;
    }
}

[BurstCompile]
public struct RopePostSolveVelocityJob : IJobParallelFor
{
    public NativeArray<float3> PointPos;
    [ReadOnly] public NativeArray<float3> PrevPos;
    public NativeArray<float3> Vel;
    public float OneOverDt;
    public float Dt;
    public float Damping; // 阻尼系数

    public void Execute(int i)
    {
        // 阻尼力
        float3 displacement = PointPos[i] - PrevPos[i];
        float dispLen = math.length(displacement);
        if (dispLen > 1e-8f && Damping > 0f)
        {
            float3 dampingCorrection = -Damping * displacement;
            PointPos[i] += dampingCorrection;
        }

        Vel[i] = (PointPos[i] - PrevPos[i]) * OneOverDt;
    }
}

[BurstCompile]
public struct RopePostSolveGhostVelJob : IJobParallelFor
{
    [ReadOnly] public NativeArray<float3> Vel;
    public NativeArray<float3> GhostVels;

    public void Execute(int i)
    {
        // i的范围是 [0, NumPoints-2)，即GhostVels的长度
        GhostVels[i] = (Vel[i] + Vel[i + 1]) / 2f;
    }
}

[BurstCompile]
public struct RopeRenderingJob : IJob
{
    [ReadOnly] public NativeArray<float3> PointPos;
    [ReadOnly] public NativeArray<float3> GhostPos;
    public NativeArray<float3> OutputVertices;
    public int Subdivision;
    public float Radius;
    public int NumPoints;

    public void Execute()
    {
        int writeIndex = 0;
        OutputVertices[writeIndex++] = PointPos[0];

        int interpolations = 10;

        for (int e = 0; e < NumPoints - 2; e++)
        {
            float3 v1 = PointPos[e];
            float3 v2 = PointPos[e + 1];
            float3 v3 = PointPos[e + 2];

            float3 vm = 0.5f * (v1 + v2);
            float3 vml = 0.5f * (v2 + v3);

            float3 d3f_raw = v2 - v1;
            float d3f_len = math.length(d3f_raw);
            if (d3f_len < 1e-8f) continue; // 跳过退化段
            float3 d3f = d3f_raw / d3f_len;

            float3 ghostDiff1 = GhostPos[e] - v1;
            float ghostDiff1Len = math.length(ghostDiff1);
            if (ghostDiff1Len < 1e-8f) continue;
            float3 crossGhost = math.cross(d3f, ghostDiff1 / ghostDiff1Len);
            float crossGhostLen = math.length(crossGhost);
            if (crossGhostLen < 1e-8f) continue;
            float3 d2f = crossGhost / crossGhostLen;
            float3 d1f = math.cross(d2f, d3f);

            float3 d3l_raw = v3 - v2;
            float d3l_len = math.length(d3l_raw);
            if (d3l_len < 1e-8f) continue;
            float3 d3l = d3l_raw / d3l_len;

            float3 ghostDiff2 = GhostPos[e + 1] - v2;
            float ghostDiff2Len = math.length(ghostDiff2);
            if (ghostDiff2Len < 1e-8f) continue;
            float3 crossGhost2 = math.cross(d3l, ghostDiff2 / ghostDiff2Len);
            float crossGhost2Len = math.length(crossGhost2);
            if (crossGhost2Len < 1e-8f) continue;
            float3 d2l = crossGhost2 / crossGhost2Len;
            float3 d1l = math.cross(d2l, d3l);

            // 构建旋转矩阵
            float3x3 De1 = new float3x3(d1f, d2f, d3f);
            float3x3 De2 = new float3x3(d1l, d2l, d3l);

            float3x3 rotation = math.mul(De2, math.transpose(De1));
            float trace = rotation.c0.x + rotation.c1.y + rotation.c2.z;
            float theta = math.acos(math.clamp((trace - 1f) / 2f, -1f, 1f));

            float3 n = float3.zero;
            if (theta > 1e-6f)
            {
                float sinTheta = math.sin(theta);
                n = new float3(
                    (rotation.c1.z - rotation.c2.y) / (2f * sinTheta),
                    (rotation.c2.x - rotation.c0.z) / (2f * sinTheta),
                    (rotation.c0.y - rotation.c1.x) / (2f * sinTheta)
                );
            }

            float le = math.length(vml - vm);
            if (le < 1e-8f) continue; // 跳过退化段
            float segmentLength = le / (interpolations - 1);

            // 插值帧
            float3 currentPosition = vm;

            for (int i = 0; i < interpolations - 1; i++)
            {
                float3x3 interpolatedFrame;
                if (i == 0)
                {
                    interpolatedFrame = De1;
                    currentPosition += d3f * segmentLength;
                }
                else
                {
                    float r = (i * segmentLength) / le;
                    float nLen = math.length(n);
                    float3 scaledTheta = nLen > 1e-8f ? theta * (n / nLen) * r : float3.zero;
                    float mag = math.length(scaledTheta);
                    float cosR = math.cos(mag);
                    float sinR = math.sin(mag);
                    float3 axis = mag > 1e-8f ? scaledTheta / mag : new float3(1, 0, 0);

                    // Rodrigues旋转公式构建矩阵
                    float3x3 rotInterp = new float3x3(
                        new float3(cosR + axis.x * axis.x * (1 - cosR),
                                   axis.y * axis.x * (1 - cosR) + axis.z * sinR,
                                   axis.z * axis.x * (1 - cosR) - axis.y * sinR),
                        new float3(axis.x * axis.y * (1 - cosR) - axis.z * sinR,
                                   cosR + axis.y * axis.y * (1 - cosR),
                                   axis.z * axis.y * (1 - cosR) + axis.x * sinR),
                        new float3(axis.x * axis.z * (1 - cosR) + axis.y * sinR,
                                   axis.y * axis.z * (1 - cosR) - axis.x * sinR,
                                   cosR + axis.z * axis.z * (1 - cosR))
                    );

                    interpolatedFrame = math.mul(rotInterp, De1);
                    float3 d3raw = interpolatedFrame.c2;
                    float d3rawLen = math.length(d3raw);
                    float3 d3 = d3rawLen > 1e-8f ? d3raw / d3rawLen : d3f;
                    currentPosition += d3 * segmentLength;
                }

                // 生成截面顶点
                for (int j = 0; j < Subdivision; j++)
                {
                    float angle = 2f * math.PI * j / Subdivision;
                    float3 frameC0 = interpolatedFrame.c0;
                    float frameC0Len = math.length(frameC0);
                    float3 frameC1 = interpolatedFrame.c1;
                    float frameC1Len = math.length(frameC1);
                    float3 vertex = currentPosition +
                        Radius * (math.cos(angle) * (frameC0Len > 1e-8f ? frameC0 / frameC0Len : new float3(1,0,0)) +
                                  math.sin(angle) * (frameC1Len > 1e-8f ? frameC1 / frameC1Len : new float3(0,1,0)));

                    if (writeIndex < OutputVertices.Length)
                        OutputVertices[writeIndex++] = vertex;
                }
            }
        }

        // 最后一个点
        if (writeIndex < OutputVertices.Length)
            OutputVertices[writeIndex] = PointPos[NumPoints - 1];
    }
}
