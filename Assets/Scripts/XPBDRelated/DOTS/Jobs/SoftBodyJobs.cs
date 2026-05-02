using Unity.Burst;
using Unity.Collections;
using Unity.Jobs;
using Unity.Mathematics;

// ============================================================================
// 软体专属 Job：四面体体积约束
// 说明：通用的 PreSolve / PostSolve / ResetLambda / DistanceConstraint /
//       AnalyticalCollision 已迁移至 XPBDCommonJobs.cs（与布料共用）。
// ============================================================================

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
