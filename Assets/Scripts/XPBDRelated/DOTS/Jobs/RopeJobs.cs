using Unity.Burst;
using Unity.Collections;
using Unity.Jobs;
using Unity.Mathematics;

// ============================================================
// Rope PreSolve Job - 重力积分 + 速度更新（可并行）
// ============================================================
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

// ============================================================
// Rope Lambda重置 Job
// ============================================================
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

// ============================================================
// Rope EdgeConstraint Job（顺序依赖，使用IJob + Burst）
// ============================================================
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

    public void Execute()
    {
        const float EPSILON = 1e-6f;
        float alpha = Stiffness / (Dt * Dt);

        for (int i = 0; i < NumPoints - 1; i++)
        {
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
                float deltaLambda = -(math.dot(p2pm, p1p0) + alpha * Lambda1[i]) / (wSum + alpha);
                PointPos[i] += math.normalize(p0p2) * deltaLambda * PointInvMass[i];
                PointPos[i + 1] += math.normalize(p2p1) * deltaLambda * PointInvMass[i + 1];
                GhostPos[i] += math.normalize(p1p0) * deltaLambda * GhostInvMass[i];
                Lambda1[i] += deltaLambda;
            }

            // === 约束3：Ghost点距离约束 ===
            wSum = 0.25f * PointInvMass[i] + 0.25f * PointInvMass[i + 1] + 1.0f * GhostInvMass[i];

            if (wSum > EPSILON)
            {
                pm = 0.5f * (PointPos[i] + PointPos[i]); // 注意：原代码就是这样写的
                p2pm = GhostPos[i] - pm;

                float p2pm_mag = math.length(p2pm);
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

// ============================================================
// Rope BendingAndTwisting Constraint Job（顺序依赖，使用IJob + Burst）
// ============================================================
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

    public void Execute()
    {
        float alpha = Stiffness / (Dt * Dt);

        for (int idx = 0; idx < NumPoints - 2; idx++)
        {
            float halfLength = RestLengths[0] / 2f;

            // 计算材料坐标系
            float3x3 a = ComputeMaterialFrame(PointPos[idx], PointPos[idx + 1], GhostPos[idx]);
            float3x3 b = ComputeMaterialFrame(PointPos[idx + 1], PointPos[idx + 2], GhostPos[idx + 1]);

            // 计算Darboux向量
            float3 darboux = ComputeDarbouxVector(a, b, halfLength);

            // 计算材料坐标系导数 dajpi[axis][point], dbjpi[axis][point]
            // axis: 0=d1, 1=d2, 2=d3; point: 0=p0, 1=p1, 2=p2(ghost)
            // 展开为独立变量以避免数组分配
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
        float3 d3 = math.normalize(p2 - p1);
        float3 d2 = math.normalize(math.cross(d3, pg - p1));
        float3 d1 = math.cross(d2, d3);
        // 列存储：c0=d1, c1=d2, c2=d3
        return new float3x3(d1, d2, d3);
    }

    static float3 ComputeDarbouxVector(float3x3 dA, float3x3 dB, float halfLength)
    {
        float factor = 1.0f + math.dot(dA.c0, dB.c0) + math.dot(dA.c1, dB.c1) + math.dot(dA.c2, dB.c2);
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
        x = 2.0f / (halfLength * x);

        // 排列组合索引
        // c=0: i=0, j=2, k=1
        // c=1: i=1, j=0, k=2
        // c=2: i=2, j=1, k=0

        // 获取dajpi和dbjpi的列（按axis索引）
        // dajpi[axis][point]: da{axis}p{point}
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

// ============================================================
// Rope PostSolve - 分两步：先更新速度+阻尼（并行），再更新GhostVels（并行，因为Vel已全部写完）
// ============================================================
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

// ============================================================
// Rope Cosserat渲染插值 Job（Burst加速）
// ============================================================
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

            float3 d3f = math.normalize(v2 - v1);
            float3 crossGhost = math.cross(d3f, math.normalize(GhostPos[e] - v1));
            float3 d2f = math.normalize(crossGhost);
            float3 d1f = math.normalize(math.cross(d2f, d3f));

            float3 d3l = math.normalize(v3 - v2);
            float3 crossGhost2 = math.cross(d3l, math.normalize(GhostPos[e + 1] - v2));
            float3 d2l = math.normalize(crossGhost2);
            float3 d1l = math.normalize(math.cross(d2l, d3l));

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
                    float3 scaledTheta = theta * math.normalize(n) * r;
                    float mag = math.length(scaledTheta);
                    float cosR = math.cos(mag);
                    float sinR = math.sin(mag);
                    float3 axis = math.normalize(scaledTheta);

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
                    float3 d3 = math.normalize(interpolatedFrame.c2);
                    currentPosition += d3 * segmentLength;
                }

                // 生成截面顶点
                for (int j = 0; j < Subdivision; j++)
                {
                    float angle = 2f * math.PI * j / Subdivision;
                    float3 vertex = currentPosition +
                        Radius * (math.cos(angle) * math.normalize(interpolatedFrame.c0) +
                                  math.sin(angle) * math.normalize(interpolatedFrame.c1));

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
