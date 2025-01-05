using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using UnityEngine;

public class EdgeConstraint : Constraint
{
    public RopeXPBDSolver solver;
    float[] length;
    const float EPSILON = 1e-6f;
    float alpha;
    private float[] lambdas0;
    private float[] lambdas1;
    private float[] lambdas2;

    public EdgeConstraint(RopeXPBDSolver solver) : base(solver)
    {
        stiff = 0.001f;
        this.solver = solver;
        this.length = solver.length;

        ResetLambda();
    }

    public override void ResetLambda()
    {
        lambdas0 = new float[solver.pointPos.Count()];
        lambdas1 = new float[solver.pointPos.Count()];
        lambdas2 = new float[solver.pointPos.Count()];
    }

    public override void SolveConstraint(float dt)
    {
        float deltaLambda;

        alpha = stiff / (dt * dt);

        for (int i = 0; i < solver.pointPos.Count() - 1; i++)
        {
            Vector3 dir = solver.pointPos[i] - solver.pointPos[i + 1];
            float len = dir.magnitude;
            float wSum = solver.pointInvMass[i] + solver.pointInvMass[i + 1];
            if (len > EPSILON && wSum > EPSILON)
            {
                //deltaLambda = -((len - length[i])) / (wSum) * edgeKs;
                deltaLambda = -((len - length[i]) + alpha * lambdas0[i]) / (wSum + alpha);

                Vector3 dP = deltaLambda * (dir / len);
                solver.pointPos[i] += dP * solver.pointInvMass[i];
                solver.pointPos[i + 1] -= dP * solver.pointInvMass[i + 1];

                lambdas0[i] += deltaLambda;
            }


            Vector3 pm = 0.5f * (solver.pointPos[i] + solver.pointPos[i + 1]);
            Vector3 p0p2 = solver.pointPos[i] - solver.ghostPos[i];
            Vector3 p2p1 = solver.ghostPos[i] - solver.pointPos[i + 1];
            Vector3 p1p0 = solver.pointPos[i + 1] - solver.pointPos[i];
            Vector3 p2pm = solver.ghostPos[i] - pm;

            wSum = solver.pointInvMass[i] * p0p2.sqrMagnitude + solver.pointInvMass[i + 1] * p2p1.sqrMagnitude + solver.ghostInvMass[i] * p1p0.sqrMagnitude;
            if (wSum > EPSILON)
            {
                deltaLambda = -(Vector3.Dot(p2pm, p1p0) + alpha * lambdas1[i]) / (wSum + alpha);
                //deltaLambda = -(Vector3.Dot(p2pm, p1p0)) / (wSum) * edgeKs;

                solver.pointPos[i] += p0p2.normalized * deltaLambda * solver.pointInvMass[i];
                solver.pointPos[i + 1] += p2p1.normalized * deltaLambda * solver.pointInvMass[i + 1];
                solver.ghostPos[i] += p1p0.normalized * deltaLambda * solver.ghostInvMass[i];

                lambdas1[i] += deltaLambda;
            }

            wSum = 0.25f * solver.pointInvMass[i] + 0.25f * solver.pointInvMass[i + 1] + 1.0f * solver.ghostInvMass[i];

            if (wSum > EPSILON)
            {
                pm = 0.5f * (solver.pointPos[i] + solver.pointPos[i]);
                p2pm = solver.ghostPos[i] - pm;

                float p2pm_mag = p2pm.magnitude;
                p2pm *= 1.0f / p2pm_mag;

                deltaLambda = -(p2pm_mag - solver.ghostDistance + alpha * lambdas2[i]) / (wSum + alpha);
                //deltaLambda = -(p2pm_mag - solver.ghostDistance) / (wSum) * edgeKs;

                solver.pointPos[i] -= 0.5f * solver.pointInvMass[i] * deltaLambda * p2pm;
                solver.pointPos[i + 1] -= 0.5f * solver.pointInvMass[i + 1] * deltaLambda * p2pm;
                solver.ghostPos[i] += 1.0f * solver.ghostInvMass[i] * deltaLambda * p2pm;

                lambdas2[i] += deltaLambda;
            }
        }
    }
}


public class BendingAndTwistingConstraint : Constraint
{
    int[][] permutation = new int[][]{
        new int[]{ 0, 2, 1 },
        new int[]{ 1, 0, 2 },
        new int[]{ 2, 1, 0 }
    };
    float bendingAndTwistingKs = 0.8f;

    public RopeXPBDSolver solver;
    Vector3[] initDarboux;
    float alpha;
    private Vector3[] lambdas;

    public BendingAndTwistingConstraint(RopeXPBDSolver solver) : base(solver)
    {
        this.solver = solver;

        initDarboux = new Vector3[solver.pointPos.Count() - 2];
        for (int i = 0; i < solver.pointPos.Count() - 2; i++)
        {
            Matrix4x4 a = ComputeMaterialFrame(solver.pointPos[i], solver.pointPos[i + 1], solver.ghostPos[i]);
            Matrix4x4 b = ComputeMaterialFrame(solver.pointPos[i + 1], solver.pointPos[i + 2], solver.ghostPos[i + 1]);

            Vector3 darboux = ComputeDarbouxVector(a, b, solver.length[0] / 2);

            initDarboux[i] = darboux;
        }

        ResetLambda();
    }

    public override void ResetLambda()
    {
        lambdas = new Vector3[solver.pointPos.Count() - 2];
    }

    public override void SolveConstraint(float dt)
    {
        alpha = stiff / (dt * dt);

        for (int i = 0; i < solver.pointPos.Count() - 2; i++)
        {
            Matrix4x4 a = ComputeMaterialFrame(solver.pointPos[i], solver.pointPos[i + 1], solver.ghostPos[i]);
            Matrix4x4 b = ComputeMaterialFrame(solver.pointPos[i + 1], solver.pointPos[i + 2], solver.ghostPos[i + 1]);

            Vector3 darboux = ComputeDarbouxVector(a, b, solver.length[0] / 2);

            Matrix4x4[][] dajpi = new Matrix4x4[3][];
            dajpi[0] = new Matrix4x4[3];
            dajpi[1] = new Matrix4x4[3];
            dajpi[2] = new Matrix4x4[3];
            ComputeMaterialFrameDerivative(solver.pointPos[i], solver.pointPos[i + 1], solver.ghostPos[i], a,
                out dajpi[0][0], out dajpi[0][1], out dajpi[0][2],
                out dajpi[1][0], out dajpi[1][1], out dajpi[1][2],
                out dajpi[2][0], out dajpi[2][1], out dajpi[2][2]);

            Matrix4x4[][] dbjpi = new Matrix4x4[3][];
            dbjpi[0] = new Matrix4x4[3];
            dbjpi[1] = new Matrix4x4[3];
            dbjpi[2] = new Matrix4x4[3];
            ComputeMaterialFrameDerivative(solver.pointPos[i + 1], solver.pointPos[i + 2], solver.ghostPos[i + 1], b,
                out dbjpi[0][0], out dbjpi[0][1], out dbjpi[0][2],
                out dbjpi[1][0], out dbjpi[1][1], out dbjpi[1][2],
                out dbjpi[2][0], out dbjpi[2][1], out dbjpi[2][2]);

            Matrix4x4[] constraint_jacobian = new Matrix4x4[5];
            ComputeDarbouxGradient(darboux, solver.length[0] / 2, a, b,
                dajpi, dbjpi,
                bendingAndTwistingKs,
                out constraint_jacobian[0],
                out constraint_jacobian[1],
                out constraint_jacobian[2],
                out constraint_jacobian[3],
                out constraint_jacobian[4]);

            Vector3 C = new Vector3((darboux[0] - initDarboux[i][0]),
                                     (darboux[1] - initDarboux[i][1]),
                                     (darboux[2] - initDarboux[i][2]));
            C *= bendingAndTwistingKs;
            Matrix4x4 factor_matrix = new Matrix4x4();
            Matrix4x4 tmp_mat = new Matrix4x4();

            float[] invMasses = { solver.pointInvMass[i], solver.pointInvMass[i + 1], solver.pointInvMass[i + 2], solver.ghostInvMass[i], solver.ghostInvMass[i + 1] };
            for (int j = 0; j < 5; ++j)
            {
                tmp_mat = Matrix4x4.Transpose(constraint_jacobian[j]) * constraint_jacobian[j];
                tmp_mat.SetColumn(0, tmp_mat.GetColumn(0) * invMasses[j]);
                tmp_mat.SetColumn(1, tmp_mat.GetColumn(1) * invMasses[j]);
                tmp_mat.SetColumn(2, tmp_mat.GetColumn(2) * invMasses[j]);

                factor_matrix = Add(factor_matrix, tmp_mat);
            }

            for (int j = 0; j < 3; j++)
            {
                factor_matrix.SetColumn(j, new Vector3(factor_matrix.GetColumn(j).x * alpha, factor_matrix.GetColumn(j).y * alpha, factor_matrix.GetColumn(j).z * alpha));
            }
            tmp_mat = Matrix4x4.Inverse(factor_matrix);

            Vector3 deltaLambda = Multiply(-(C + alpha * lambdas[i]), tmp_mat);

            lambdas[i] += deltaLambda;

            Vector3[] dp = new Vector3[5];

            for (int j = 0; j < 5; ++j)
            {
                constraint_jacobian[j].SetColumn(0, constraint_jacobian[j].GetColumn(0) * invMasses[j]);
                constraint_jacobian[j].SetColumn(1, constraint_jacobian[j].GetColumn(1) * invMasses[j]);
                constraint_jacobian[j].SetColumn(2, constraint_jacobian[j].GetColumn(2) * invMasses[j]);
                //dp[j] = (constraint_jacobian[j]) * (tmp_mat * C);
                dp[j] = (constraint_jacobian[j]) * deltaLambda;
            }

            solver.pointPos[i] += dp[0];
            solver.pointPos[i + 1] += dp[1];
            solver.pointPos[i + 2] += dp[2];
            solver.ghostPos[i] += dp[3];
            solver.ghostPos[i + 1] += dp[4];
        }
    }

    Matrix4x4 ComputeMaterialFrame(Vector3 p1, Vector3 p2, Vector3 pg)
    {
        Matrix4x4 frame = new Matrix4x4();
        frame.SetColumn(2, p2 - p1);
        frame.GetColumn(2).Normalize();

        frame.SetColumn(1, Vector3.Cross(frame.GetColumn(2), pg - p1));
        frame.GetColumn(1).Normalize();

        frame.SetColumn(0, Vector3.Cross(frame.GetColumn(1), frame.GetColumn(2)));
        //	frame.col(0).normalize();
        return frame;
    }

    Vector3 ComputeDarbouxVector(Matrix4x4 dA, Matrix4x4 dB, float halfLength)
    {
        Vector3 darboux = new Vector3();
        float factor = 1.0f + Vector3.Dot(dA.GetColumn(0), dB.GetColumn(0)) + Vector3.Dot(dA.GetColumn(1), dB.GetColumn(1)) + Vector3.Dot(dA.GetColumn(2), dB.GetColumn(2));

        factor = 2.0f / (halfLength * factor);

        for (int c = 0; c < 3; ++c)
        {
            int i = permutation[c][0];
            int j = permutation[c][1];
            int k = permutation[c][2];

            darboux[i] = Vector3.Dot(dA.GetColumn(j), dB.GetColumn(k)) - Vector3.Dot(dA.GetColumn(k), dB.GetColumn(j));
        }

        darboux *= factor;
        return darboux;
    }

    void ComputeMaterialFrameDerivative(Vector3 p0, Vector3 p1, Vector3 p2, Matrix4x4 d,
    out Matrix4x4 d1p0, out Matrix4x4 d1p1, out Matrix4x4 d1p2,
    out Matrix4x4 d2p0, out Matrix4x4 d2p1, out Matrix4x4 d2p2,
    out Matrix4x4 d3p0, out Matrix4x4 d3p1, out Matrix4x4 d3p2)
    {
        // d3pi
        Vector3 p01 = p1 - p0;
        float length_p01 = p01.magnitude;
        d3p0 = new Matrix4x4();
        d3p2 = new Matrix4x4();
        d3p1 = new Matrix4x4();

        d3p0.SetColumn(0, d.GetColumn(2)[0] * d.GetColumn(2));
        d3p0.SetColumn(1, d.GetColumn(2)[1] * d.GetColumn(2));
        d3p0.SetColumn(2, d.GetColumn(2)[2] * d.GetColumn(2));

        d3p0.SetColumn(0, new Vector3(d3p0.GetColumn(0).x - 1, d3p0.GetColumn(0).y, d3p0.GetColumn(0).z));
        d3p0.SetColumn(1, new Vector3(d3p0.GetColumn(1).x, d3p0.GetColumn(0).y - 1, d3p0.GetColumn(0).z));
        d3p0.SetColumn(2, new Vector3(d3p0.GetColumn(2).x, d3p0.GetColumn(0).y, d3p0.GetColumn(0).z - 1));

        d3p0.SetColumn(0, d3p0.GetColumn(0) * (1.0f / length_p01));
        d3p0.SetColumn(1, d3p0.GetColumn(1) * (1.0f / length_p01));
        d3p0.SetColumn(2, d3p0.GetColumn(2) * (1.0f / length_p01));

        d3p1.SetColumn(0, -d3p0.GetColumn(0));
        d3p1.SetColumn(1, -d3p0.GetColumn(1));
        d3p1.SetColumn(2, -d3p0.GetColumn(2));

        d3p2.SetColumn(0, new Vector4());
        d3p2.SetColumn(1, new Vector4());
        d3p2.SetColumn(2, new Vector4());

        ////>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        //// d2pi
        Vector3 p02 = p2 - p0;
        Vector3 p01_cross_p02 = Vector3.Cross(p01, p02);

        float length_cross = p01_cross_p02.magnitude;

        Matrix4x4 mat = new Matrix4x4();
        mat.SetColumn(0, d.GetColumn(1)[0] * d.GetColumn(1));
        mat.SetColumn(1, d.GetColumn(1)[1] * d.GetColumn(1));
        mat.SetColumn(2, d.GetColumn(1)[2] * d.GetColumn(1));

        mat.SetColumn(0, new Vector3(mat.GetColumn(0).x - 1, mat.GetColumn(0).y, mat.GetColumn(0).z));
        mat.SetColumn(1, new Vector3(mat.GetColumn(1).x, mat.GetColumn(0).y - 1, mat.GetColumn(0).z));
        mat.SetColumn(2, new Vector3(mat.GetColumn(2).x, mat.GetColumn(0).y, mat.GetColumn(0).z - 1));

        mat.SetColumn(0, mat.GetColumn(0) * (-1.0f / length_cross));
        mat.SetColumn(1, mat.GetColumn(1) * (-1.0f / length_cross));
        mat.SetColumn(2, mat.GetColumn(2) * (-1.0f / length_cross));

        Matrix4x4 product_matrix = Cross(p2 - p1);
        d2p0 = mat * product_matrix;

        product_matrix = Cross(p0 - p2);
        d2p1 = mat * product_matrix;

        product_matrix = Cross(p1 - p0);
        d2p2 = mat * product_matrix;

        ////>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        //// d1pi
        Matrix4x4 product_mat_d3;
        Matrix4x4 product_mat_d2;
        Matrix4x4 m1, m2;

        product_mat_d3 = Cross(d.GetColumn(2));
        product_mat_d2 = Cross(d.GetColumn(1));

        //dj::MulMatrix3x3<real>(&product_mat_d3[0][0], &d2p0[0][0], &m1[0][0]);
        //dj::MulMatrix3x3<real>(&product_mat_d2[0][0], &d3p0[0][0], &m2[0][0]);
        d1p0 = Minus(product_mat_d2 * d3p0, product_mat_d3 * d2p0);

        //dj::MulMatrix3x3<real>(&product_mat_d3[0][0], &d2p1[0][0], &m1[0][0]);
        //dj::MulMatrix3x3<real>(&product_mat_d2[0][0], &d3p1[0][0], &m2[0][0]);
        d1p1 = Minus(product_mat_d2 * d3p1, product_mat_d3 * d2p1);

        /*dj::MulMatrix3x3<real>(&product_mat_d3[0][0], &d2p2[0][0], &d1p2[0][0]);*/
        d1p2 = product_mat_d3 * d2p2;
        d1p2.SetColumn(0, d1p2.GetColumn(0) * (-1.0f));
        d1p2.SetColumn(1, d1p2.GetColumn(1) * (-1.0f));
        d1p2.SetColumn(2, d1p2.GetColumn(2) * (-1.0f));
    }
    //0,     -v.Z,  v.Y,
    //v.Z,    0,    -v.X,
    //-v.Y,  v.X,      0
    Matrix4x4 Cross(Vector3 mat)
    {
        Matrix4x4 newMat = new Matrix4x4();

        newMat.SetColumn(0, new Vector3(0, mat.z, -mat.y));
        newMat.SetColumn(1, new Vector3(-mat.z, 0, mat.x));
        newMat.SetColumn(2, new Vector3(mat.y, -mat.x, 0));

        return newMat;
    }

    Matrix4x4 Minus(Matrix4x4 mat1, Matrix4x4 mat2)
    {
        Matrix4x4 newMat = new Matrix4x4();

        newMat.SetColumn(0, mat1.GetColumn(0) - mat2.GetColumn(0));
        newMat.SetColumn(1, mat1.GetColumn(1) - mat2.GetColumn(1));
        newMat.SetColumn(2, mat1.GetColumn(2) - mat2.GetColumn(2));

        return newMat;
    }

    Matrix4x4 Add(Matrix4x4 mat1, Matrix4x4 mat2)
    {
        Matrix4x4 newMat = new Matrix4x4();

        newMat.SetColumn(0, mat1.GetColumn(0) + mat2.GetColumn(0));
        newMat.SetColumn(1, mat1.GetColumn(1) + mat2.GetColumn(1));
        newMat.SetColumn(2, mat1.GetColumn(2) + mat2.GetColumn(2));

        return newMat;
    }
    //v.X* m.M11 + v.Y* m.M21 + v.Z* m.M31,
    //v.X* m.M12 + v.Y* m.M22 + v.Z* m.M32,
    //v.X* m.M13 + v.Y* m.M23 + v.Z* m.M33
    Vector3 Multiply(Vector3 vec, Matrix4x4 mat)
    {
        float[] num = new float[3];
        for (int i = 0; i < 3; i++)
        {
            num[i] += vec.x * mat.GetColumn(i).x + vec.y * mat.GetColumn(i).y + vec.z * mat.GetColumn(i).z;
        }
        return new Vector3(num[0], num[1], num[2]);
    }

    void ComputeDarbouxGradient(
    Vector3 darboux_vector, float length,
    Matrix4x4 da, Matrix4x4 db,
    Matrix4x4[][] dajpi, Matrix4x4[][] dbjpi,
    float bendAndTwistKs,
    out Matrix4x4 omega_pa, out Matrix4x4 omega_pb, out Matrix4x4 omega_pc, out Matrix4x4 omega_pd, out Matrix4x4 omega_pe
    )
    {
        omega_pa = new Matrix4x4();
        omega_pb = new Matrix4x4();
        omega_pc = new Matrix4x4();
        omega_pd = new Matrix4x4();
        omega_pe = new Matrix4x4();

        //float x = 1.0f + da[0] * db[0] + da[1] * db[1] + da[2] * db[2];
        float x = 1.0f + Vector3.Dot(da.GetColumn(0), db.GetColumn(0)) + Vector3.Dot(da.GetColumn(1), db.GetColumn(1)) + Vector3.Dot(da.GetColumn(2), db.GetColumn(2));
        x = 2.0f / (length * x);

        for (int c = 0; c < 3; ++c)
        {
            int i = permutation[c][0];
            int j = permutation[c][1];
            int k = permutation[c][2];
            // pa
            Vector3 term1 = new Vector3(0, 0, 0);
            Vector3 term2 = new Vector3(0, 0, 0);
            Vector3 tmp = new Vector3(0, 0, 0);

            {

                //DOUBLE CHECK !!!
                term1 = Matrix4x4.Transpose(dajpi[j][0]) * db.GetColumn(k);
                tmp = Matrix4x4.Transpose(dajpi[k][0]) * db.GetColumn(j);
                term1 = term1 - tmp;
                // second term
                for (int n = 0; n < 3; ++n)
                {
                    tmp = Matrix4x4.Transpose(dajpi[n][0]) * db.GetColumn(n);
                    term2 = term2 + tmp;
                }
                omega_pa = new Matrix4x4();
                omega_pa.SetColumn(i, (term1) - (0.5f * darboux_vector[i] * length) * (term2));
                omega_pa.SetColumn(i, omega_pa.GetColumn(i) * x * bendAndTwistKs);
            }
            // pb
            {
                term1 = new(0, 0, 0);
                term2 = new(0, 0, 0);
                tmp = new(0, 0, 0);
                // first term
                //dj::MulVecMatrix3x3<real>(db[k](), (real(*)[3]) &dajpi[j][1], term1());
                //dj::MulVecMatrix3x3<real>(db[j](), (real(*)[3]) &dajpi[k][1], tmp());
                term1 = Matrix4x4.Transpose(dajpi[j][1]) * db.GetColumn(k);
                tmp = Matrix4x4.Transpose(dajpi[k][1]) * db.GetColumn(j);
                term1 = term1 - tmp;
                // third term
                //dj::MulVecMatrix3x3<real>(da[k](), (real(*)[3]) &dbjpi[j][0], tmp());
                tmp = Matrix4x4.Transpose(dbjpi[j][0]) * da.GetColumn(k);
                term1 = term1 - tmp;

                //dj::MulVecMatrix3x3<real>(da[j](), (real(*)[3]) &dbjpi[k][0], tmp());
                tmp = Matrix4x4.Transpose(dbjpi[k][0]) * da.GetColumn(j);
                term1 = term1 + tmp;

                // second term
                for (int n = 0; n < 3; ++n)
                {
                    //dj::MulVecMatrix3x3<real>(db[n](), (real(*)[3]) &dajpi[n][1], tmp());
                    tmp = Matrix4x4.Transpose(dajpi[n][1]) * db.GetColumn(n);
                    term2 = term2 + tmp;

                    //dj::MulVecMatrix3x3<real>(da[n](), (real(*)[3]) &dbjpi[n][0], tmp());
                    tmp = Matrix4x4.Transpose(dbjpi[n][0]) * da.GetColumn(n);
                    term2 = term2 + tmp;
                }

                omega_pb.SetColumn(i, (term1) - (0.5f * darboux_vector[i] * length) * (term2));
                omega_pb.SetColumn(i, omega_pb.GetColumn(i) * x * bendAndTwistKs);
            }
            // pc
            {
                term1 = new(0, 0, 0);
                term2 = new(0, 0, 0);
                tmp = new(0, 0, 0);

                // first term
                //dj::MulVecMatrix3x3<real>(da[k](), (real(*)[3]) &dbjpi[j][1], term1());
                //dj::MulVecMatrix3x3<real>(da[j](), (real(*)[3]) &dbjpi[k][1], tmp());
                term1 = Matrix4x4.Transpose(dbjpi[j][1]) * da.GetColumn(k);
                tmp = Matrix4x4.Transpose(dbjpi[k][1]) * da.GetColumn(j);
                term1 = term1 - tmp;

                // second term
                for (int n = 0; n < 3; ++n)
                {
                    //dj::MulVecMatrix3x3<real>(da[n](), (real(*)[3]) &dbjpi[n][1], tmp());
                    tmp = Matrix4x4.Transpose(dbjpi[n][1]) * da.GetColumn(n);
                    term2 = term2 + tmp;
                }
                omega_pc = new Matrix4x4();
                omega_pc.SetColumn(i, (term1) + (0.5f * darboux_vector[i] * length) * (term2));
                omega_pc.SetColumn(i, omega_pc.GetColumn(i) * -x * bendAndTwistKs);
            }
            // pd
            {
                term1 = new(0, 0, 0);
                term2 = new(0, 0, 0);
                tmp = new(0, 0, 0);
                // first term
                //dj::MulVecMatrix3x3<real>(db[k](), (real(*)[3]) &dajpi[j][2], term1());
                //dj::MulVecMatrix3x3<real>(db[j](), (real(*)[3]) &dajpi[k][2], tmp());
                term1 = Matrix4x4.Transpose(dajpi[j][2]) * db.GetColumn(k);
                tmp = Matrix4x4.Transpose(dajpi[k][2]) * db.GetColumn(j);
                term1 = term1 - tmp;
                // second term
                for (int n = 0; n < 3; ++n)
                {
                    //dj::MulVecMatrix3x3<real>(db[n](), (real(*)[3]) &dajpi[n][2], tmp());
                    tmp = Matrix4x4.Transpose(dajpi[n][2]) * db.GetColumn(n);
                    term2 = term2 + tmp;
                }
                omega_pd = new Matrix4x4();
                omega_pd.SetColumn(i, (term1) - (0.5f * darboux_vector[i] * length) * (term2));
                omega_pd.SetColumn(i, omega_pd.GetColumn(i) * x * bendAndTwistKs);
            }
            // pe
            {
                term1 = new(0, 0, 0);
                term2 = new(0, 0, 0);
                tmp = new(0, 0, 0);
                // first term

                //dj::MulVecMatrix3x3<real>(da[k](), (real(*)[3]) &dbjpi[j][2], term1());
                //dj::MulVecMatrix3x3<real>(da[j](), (real(*)[3]) &dbjpi[k][2], tmp());
                term1 = Matrix4x4.Transpose(dbjpi[j][2]) * da.GetColumn(k);
                tmp = Matrix4x4.Transpose(dbjpi[k][2]) * da.GetColumn(j);
                term1 -= tmp;

                // second term
                for (int n = 0; n < 3; ++n)
                {   //WARNING!!! &dbjpi[n][2][0][0] ???
                    //dj::MulVecMatrix3x3<real>(da[n](), (real(*)[3]) &dbjpi[n][2][0][0], tmp());
                    tmp = Matrix4x4.Transpose(dbjpi[n][2]) * da.GetColumn(n);
                    term2 += tmp;
                }

                omega_pe = new Matrix4x4();
                omega_pe.SetColumn(i, (term1) + (0.5f * darboux_vector[i] * length) * (term2));
                omega_pe.SetColumn(i, omega_pe.GetColumn(i) * -x * bendAndTwistKs);
            }
        }
    }
}


//public class BendTwistConstraint : Constraint
//{
//    RopeXPBDSolver solver;
//    public Quaternion[] restQs;
//    public BendTwistConstraint(XPBDSolver solver) : base(solver)
//    {
//        this.solver = solver as RopeXPBDSolver;
//        restQs = new Quaternion[this.solver.pointPos.Count() - 1];


//        for (int i = 0; i < this.solver.pointQ.Count() - 1; i++)
//        {
//            restQs[i] = CalculateDarbQ(i);
//        }
//    }


//    public Quaternion CalculateDarbQ(int i)
//    {
//        Quaternion restDarbouxVector = new Quaternion(solver.pointQ[i].x, solver.pointQ[i].y, solver.pointQ[i].z, -solver.pointQ[i].w) * solver.pointQ[i + 1];
//        Quaternion omega_plus, omega_minus;

//        omega_plus = new Quaternion(restDarbouxVector.x + 1, restDarbouxVector.y, restDarbouxVector.z, restDarbouxVector.w);
//        omega_minus = new Quaternion(restDarbouxVector.x - 1, restDarbouxVector.y, restDarbouxVector.z, restDarbouxVector.w);
//        if (SqrMagnitude(omega_minus) > SqrMagnitude(omega_plus))
//            restDarbouxVector = new Quaternion(-restDarbouxVector.x, -restDarbouxVector.y, -restDarbouxVector.z, -restDarbouxVector.w);
//        return restDarbouxVector;
//    }

//    public float SqrMagnitude(Quaternion quaternion)
//    {
//        return quaternion.x * quaternion.x + quaternion.y * quaternion.y + quaternion.z * quaternion.z + quaternion.w * quaternion.w;
//    }

//    public override void SolveConstraint(float dt)
//    {
//        for (int i = 0; i < solver.pointQ.Count() - 1; ++i)
//        {
//            Quaternion omega = new Quaternion(solver.pointQ[i].x, solver.pointQ[i].y, solver.pointQ[i].z, -solver.pointQ[i].w) * solver.pointQ[i + 1];

//            Quaternion omega_plus;
//            omega_plus = new Quaternion(restQs[i].x + omega.x, restQs[i].y + omega.y, restQs[i].z + omega.z, restQs[i].w + omega.w);
//            omega = new Quaternion(-restQs[i].x + omega.x, -restQs[i].y + omega.y, -restQs[i].z + omega.z, -restQs[i].w + omega.w);
//            if(SqrMagnitude(omega) > SqrMagnitude(omega_plus))
//            {
//                omega = omega_plus;
//            }

//            for (int j = 0; j < 3; j++) omega[j] *= stiff / (solver.invMass[i] + solver.invMass[i + 1] + (1.0e-6f));
//            omega.w = 0.0f;    //discrete Darboux vector does not have vanishing scalar part

//            solver.pointQ[i] *= new Quaternion(omega.x * solver.invMass[i], omega.y * solver.invMass[i], omega.z * solver.invMass[i], omega.w * solver.invMass[i]);
//            solver.pointQ[i + 1] *= new Quaternion(omega.x * solver.invMass[i + 1], omega.y * solver.invMass[i + 1], omega.z * solver.invMass[i + 1], omega.w * solver.invMass[i + 1]);

//            solver.pointQ[i].Normalize();
//            solver.pointQ[i + 1].Normalize();
//        }
//    }

//    public override void ResetLambda()
//    {

//    }
//}



//public class RopeCollisionConstraint : Constraint
//{
//    float restDistance = 0.1f;

//    float[] lambdas;
//    public RopeCollisionConstraint(RopeXPBDSolver solver) : base(solver)
//    {
//        stiff = 0f;
//        lambdas = new float[mySolver.numParticles];

//    }

//    public override void ResetLambda()
//    {
//        lambdas = new float[mySolver.numParticles];
//    }


//    public override void SolveConstraint(float dt)
//    {
//        var solver = (RopeXPBDSolver)mySolver;
//        if (solver.collisions.Count <= 0)
//        {
//            return;
//        }
//        float alpha = stiff / (Mathf.Pow(dt, 2));


//        foreach (var pair in solver.collisions)
//        {
//            int index = pair.Key;
//            Vector3 distance = pair.Value;
//            //Debug.Log(distance);

//            float l = distance.magnitude;
//            float l_rest = restDistance;

//            float C = l - l_rest;
//            if (C > 0.01)
//            {
//                continue;
//            }
//            //Debug.Log("wawa" +C.ToString());

//            //(xo-x1) * (1/|x0-x1|) = gradC
//            Vector3 gradC = distance.normalized;

//            float wTot = 1;

//            //lambda because |grad_Cn|^2 = 1 because if we move a particle 1 unit, the distance between the particles also grows with 1 unit, and w = w0 + w1
//            float deltalambda = (C) / (wTot);
//            lambdas[index] += deltalambda;
//            //Move the vertices x = x + deltaX where deltaX = lambda * w * gradC
//            //Debug.Log($"{index}增加的距离为{deltalambda}");
//            solver.pointPos[index] += deltalambda * gradC;
//            //solver.collisions[index].Value += (deltalambda * gradC));
//        }
//        solver.collisions.Clear();
//    }
//}

public class DoubleDistanceConstraint : Constraint
{
    public RopeXPBDSolver solver;
    float[] length;
    const float EPSILON = 1e-6f;
    float alpha;
    private float[] lambdas0;

    public DoubleDistanceConstraint(RopeXPBDSolver solver) : base(solver)
    {
        stiff = 0.1f;
        this.solver = solver;
        this.length = solver.length;

        ResetLambda();
    }

    public override void ResetLambda()
    {
        lambdas0 = new float[2];
    }

    public override void SolveConstraint(float dt)
    {
        float deltaLambda;

        alpha = stiff / (dt * dt);

        int control = solver.ctrlIndex;

        Vector3 dirLeft = solver.pointPos[control] - solver.pointPos[0];
        Vector3 dirRight = solver.pointPos[control] - solver.pointPos[solver.pointPos.Length - 1];
        float lenl = dirLeft.magnitude;
        float lenr = dirRight.magnitude;

        if (lenl > EPSILON)
        {
            float oriLen = 0f;
            for (int i = 0; i < control - 1; i++)
            {
                oriLen += length[i];
            }
            deltaLambda = -((lenl - (oriLen) * 0.9f) + alpha * lambdas0[0]) / (1 + alpha);
            lambdas0[0] += deltaLambda;

            solver.pointPos[0] -= deltaLambda * dirLeft.normalized;
        }

        if (lenr > EPSILON)
        {
            float oriLen = 0f;
            for (int i = control; i < solver.pointPos.Length - 1; i++)
            {
                oriLen += length[i];
            }
            deltaLambda = -((lenl - (oriLen) * 0.9f) + alpha * lambdas0[1]) / (1 + alpha);
            lambdas0[1] += deltaLambda;

            solver.pointPos[1] -= deltaLambda * dirRight.normalized;
        }
    }
}