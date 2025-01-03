using System.Collections.Generic;
using System.Linq;
using TMPro;
using Unity.Mathematics;
using Unity.VisualScripting;
using UnityEngine;
using static UnityEditor.PlayerSettings;
using static UnityEditor.ShaderGraph.Internal.KeywordDependentCollection;
using static ClothXPBDSolver;

public class DistanceConstraint : Constraint
{
    int[] pointAs;
    int[] pointBs;
    float[] distance;

    float[] lambdas;
    public DistanceConstraint(ClothXPBDSolver solver) : base(solver)
    {
        stiff = 0f;
        Dictionary<KeyValuePair<int, int>, float> distances = new Dictionary<KeyValuePair<int, int>, float>();

        for (int i = 0; i < mySolver.triangles.Count(); i++)
        {
            List<int> points = mySolver.triangles[i].pointsId.ToList();
            points.Sort();
            var pair = points.SelectMany((first, index) => points.Skip(index + 1), (first, second) => new { first, second });
            foreach (var point in pair)
            {
                if (!distances.ContainsKey(new(point.first, point.second)))
                {
                    distances.Add(new(point.first, point.second), (mySolver.pos[point.first] - mySolver.pos[point.second]).magnitude);
                }
            }
        }

        pointAs = new int[distances.Count];
        pointBs = new int[distances.Count];
        distance = new float[distances.Count];
        int j = 0;
        foreach (var pair in distances)
        {
            pointAs[j] = pair.Key.Key;
            pointBs[j] = pair.Key.Value;
            distance[j] = pair.Value;
            j++;
        }

        lambdas = new float[distances.Count];
    }

    public override void ResetLambda()
    {
        lambdas = new float[distance.Length];
    }


    public override void SolveConstraint(float dt)
    {
        float alpha = stiff / (math.pow(dt, 2));
        int i = 0;
        foreach (var pair in distance)
        {
            //2 vertices per edge in the data structure, so multiply by 2 to get the correct vertex index
            int id0 = pointAs[i];
            int id1 = pointBs[i];

            float w0 = mySolver.invMass[id0];
            float w1 = mySolver.invMass[id1];

            //The current length of the edge l

            //x0-x1
            //The result is stored in grads array
            Vector3 id0_minus_id1 = mySolver.pos[id0] - mySolver.pos[id1];

            //sqrMargnitude(x0-x1)
            float l = Vector3.Magnitude(id0_minus_id1);

            //If they are at the same pos we get a divisio by 0 later so ignore
            if (l == 0f)
            {
                continue;
            }

            //(xo-x1) * (1/|x0-x1|) = gradC
            Vector3 gradC = id0_minus_id1 / l;

            float l_rest = distance[i];

            float C = l - l_rest;
            float wTot = (w1 + w0) * gradC.magnitude * gradC.magnitude;

            //lambda because |grad_Cn|^2 = 1 because if we move a particle 1 unit, the distance between the particles also grows with 1 unit, and w = w0 + w1
            float deltalambda = -(C + lambdas[i] * alpha) / (wTot + alpha);
            lambdas[i++] += deltalambda;
            //Move the vertices x = x + deltaX where deltaX = lambda * w * gradC
            mySolver.pos[id0] += deltalambda * w0 * gradC;
            mySolver.pos[id1] += -deltalambda * w1 * gradC;
        }
    }
}

public class SectionDistanceConstraint : Constraint
{
    int[] pointAs;
    int[] pointBs;
    float[] distance;

    float[] lambdas;
    public SectionDistanceConstraint(ClothXPBDSolver solver) : base(solver)
    {
        stiff = 0.1f;
        Dictionary<KeyValuePair<int, int>, float> distances = new Dictionary<KeyValuePair<int, int>, float>();

        for (int i = 0; i < ((ClothXPBDSolver)mySolver).sections.Count(); i++)
        {
            for (int k = 0; k < solver.subdivision - 1; k++)
            {
                distances.Add(new KeyValuePair<int, int>(solver.sections[i][k], solver.sections[i][k + 1]), (solver.pos[solver.sections[i][k]] - solver.pos[solver.sections[i][k + 1]]).magnitude);
                if (k != 0)
                {
                    distances.Add(new KeyValuePair<int, int>(solver.sections[i][0], solver.sections[i][k + 1]), (solver.pos[solver.sections[i][0]] - solver.pos[solver.sections[i][k + 1]]).magnitude);
                }
            }
            //distances.Add(new KeyValuePair<int, int>(solver.sections[i][0], solver.sections[i][solver.subdivision - 1]), (solver.pos[solver.sections[i][0]] - solver.pos[solver.sections[i][solver.subdivision - 1]]).magnitude);
        }
        pointAs = new int[distances.Count];
        pointBs = new int[distances.Count];
        distance = new float[distances.Count];
        int j = 0;
        foreach (var pair in distances)
        {
            pointAs[j] = pair.Key.Key;
            pointBs[j] = pair.Key.Value;
            distance[j] = pair.Value;
            //Debug.Log($"{pointAs[j]}, {pointBs[j]}");
            j++;
        }

        lambdas = new float[distances.Count];
    }

    public override void ResetLambda()
    {
        lambdas = new float[distance.Length];
    }


    public override void SolveConstraint(float dt)
    {
        float alpha = stiff / (math.pow(dt, 2));
        int i = 0;
        foreach (var pair in distance)
        {
            //2 vertices per edge in the data structure, so multiply by 2 to get the correct vertex index
            int id0 = pointAs[i];
            int id1 = pointBs[i];

            float w0 = mySolver.invMass[id0];
            float w1 = mySolver.invMass[id1];

            //The current length of the edge l

            //x0-x1
            //The result is stored in grads array
            Vector3 id0_minus_id1 = mySolver.pos[id0] - mySolver.pos[id1];

            //sqrMargnitude(x0-x1)
            float l = Vector3.Magnitude(id0_minus_id1);

            //If they are at the same pos we get a divisio by 0 later so ignore
            if (l == 0f)
            {
                continue;
            }

            //(xo-x1) * (1/|x0-x1|) = gradC
            Vector3 gradC = id0_minus_id1 / l;

            float l_rest = distance[i];

            float C = l - l_rest;
            float wTot = (w1 + w0) * gradC.magnitude * gradC.magnitude;

            //lambda because |grad_Cn|^2 = 1 because if we move a particle 1 unit, the distance between the particles also grows with 1 unit, and w = w0 + w1
            float deltalambda = -(C + lambdas[i] * alpha) / (wTot + alpha);
            lambdas[i++] += deltalambda;
            //Move the vertices x = x + deltaX where deltaX = lambda * w * gradC
            mySolver.pos[id0] += deltalambda * w0 * gradC;
            mySolver.pos[id1] += -deltalambda * w1 * gradC;
        }
    }
}

public class RoundConstraint : Constraint
{
    float[] lambdas;
    float[] diameterLambdas;
    float stiff = 0.8f;
    float diameterStiff = 0.1f;
    float idealR = 0f;

    // 存储直径点对
    List<(int, int)> diameterPairs;

    public RoundConstraint(ClothXPBDSolver solver) : base(solver)
    {
        lambdas = new float[solver.numParticles];
        diameterLambdas = new float[solver.numParticles];
        lambdas = Enumerable.Repeat(0f, lambdas.Length).ToArray();
        diameterLambdas = Enumerable.Repeat(0f, lambdas.Length).ToArray();
        diameterPairs = new List<(int, int)>();

        var section = ((ClothXPBDSolver)mySolver).sections[0];
        Vector3 center = Vector3.zero;
        foreach (int pointId in section)
        {
            center += mySolver.pos[pointId];
        }
        center /= section.Length;
        idealR = (mySolver.pos[section[0]] - center).magnitude;

        // 初始化直径点对
        InitializeDiameterPairs(section);
    }

    private void InitializeDiameterPairs(int[] section)
    {
        // 为前半部分的每个点找到最适合的对应点
        for (int i = 0; i < section.Length / 2; i++)
        {
            int pointId1 = section[i];
            float maxDist = 0f;
            int bestMatch = -1;

            // 在后半部分寻找最远的点
            for (int j = section.Length / 2; j < section.Length; j++)
            {
                int pointId2 = section[j];
                float dist = (mySolver.pos[pointId2] - mySolver.pos[pointId1]).sqrMagnitude;
                if (dist > maxDist)
                {
                    maxDist = dist;
                    bestMatch = pointId2;
                }
            }

            if (bestMatch != -1)
            {
                diameterPairs.Add((pointId1, bestMatch));
            }
        }
    }

    public override void ResetLambda()
    {
        foreach (var triangle in mySolver.triangles)
        {
            triangle.lambda = 0;
        }
        for (int i = 0; i < lambdas.Length; i++)
        {
            lambdas[i] = 0;
            diameterLambdas[i] = 0;
        }
    }

    public override void SolveConstraint(float dt)
    {
        float alpha = stiff / (dt * dt);
        float diameterAlpha = diameterStiff / (dt * dt);

        foreach (var section in ((ClothXPBDSolver)mySolver).sections)
        {
            // 计算圆心
            Vector3 center = Vector3.zero;
            foreach (int pointId in section)
            {
                center += mySolver.pos[pointId];
            }
            center /= section.Length;

            // 第一步：圆形约束
            foreach (int pointId in section)
            {
                if (mySolver.invMass[pointId] == 0f) continue;

                Vector3 pointToCenter = mySolver.pos[pointId] - center;
                float currentDistance = pointToCenter.magnitude;

                if (currentDistance < 1e-6f) continue;

                Vector3 gradC = pointToCenter.normalized;
                float C = currentDistance - idealR;
                float w = mySolver.invMass[pointId];
                float wTot = w * gradC.sqrMagnitude;

                float deltalambda = -(C + lambdas[pointId] * alpha) / (wTot + alpha);
                deltalambda = Mathf.Clamp(deltalambda, -1f, 1f);
                lambdas[pointId] += deltalambda;

                mySolver.pos[pointId] += deltalambda * w * gradC;
            }

            // 第二步：直径约束
            foreach (var (pointId1, pointId2) in diameterPairs)
            {
                float w1 = mySolver.invMass[pointId1];
                float w2 = mySolver.invMass[pointId2];

                //if (w1 == 0f && w2 == 0f) continue;

                Vector3 dir = mySolver.pos[pointId2] - mySolver.pos[pointId1];
                float currentDiameter = dir.magnitude;

                //if (currentDiameter < 1e-6f) continue;

                dir = dir.normalized;
                float targetDiameter = idealR * 2f;
                float C = currentDiameter - targetDiameter;

                float wTot = (w1 + w2) * dir.sqrMagnitude;
                //if (wTot < 1e-6f) continue;

                float deltalambda = -(C + diameterLambdas[pointId1] * diameterAlpha) / (wTot + diameterAlpha);
                deltalambda = Mathf.Clamp(deltalambda, -0.5f, 0.5f);
                diameterLambdas[pointId1] += deltalambda;

                if (w1 > 0f)
                    mySolver.pos[pointId1] -= deltalambda * w1 * dir;
                if (w2 > 0f)
                    mySolver.pos[pointId2] += deltalambda * w2 * dir;
            }
        }
    }
}

public class VolumeConstraint : Constraint
{
    public VolumeConstraint(XPBDSolver solver) : base(solver)
    {
        for (int i = 0; i < solver.triangles.Length; i++)
        {
            solver.triangles[i].centroid = CalculateCentroid(solver.triangles[i]);
            solver.triangles[i].initialMatrix = CalculateShapeMatrix(solver.triangles[i]);
        }
    }

    public override void ResetLambda()
    {
        foreach (var triangle in mySolver.triangles)
        {
            triangle.lambda = 0;
        }
    }

    public Vector3 CalculateCentroid(Triangle triangle)
    {
        int index1 = triangle.pointsId[0];
        int index2 = triangle.pointsId[1];
        int index3 = triangle.pointsId[2];

        float mass1 = 1;//mySolver.invMass[index1] != 0 ? 1 / mySolver.invMass[index1] : 0;
        float mass2 = 1;// mySolver.invMass[index2] != 0 ? 1 / mySolver.invMass[index2] : 0;
        float mass3 = 1;// mySolver.invMass[index3] != 0 ? 1 / mySolver.invMass[index3] : 0;

        return (mySolver.pos[index1] * 100 * mass1 + mySolver.pos[index2] * 100 * mass2 + mySolver.pos[index3] * 100 * mass3) / (mass1 + mass2 + mass3);
    }

    // 计算形状矩阵
    public Matrix4x4 CalculateShapeMatrix(Triangle triangle)
    {
        Vector3 centroid = CalculateCentroid(triangle);
        // 三角形下标 还要点的坐标
        int index1 = triangle.pointsId[0];
        int index2 = triangle.pointsId[1];
        int index3 = triangle.pointsId[2];

        Vector3 x1_x0 = mySolver.pos[index1] * 100 - centroid;
        Vector3 x2_x0 = mySolver.pos[index2] * 100 - centroid;
        Vector3 x3_x0 = mySolver.pos[index3] * 100 - centroid;

        Matrix4x4 shapeMatrix = new Matrix4x4();
        shapeMatrix.SetRow(0, new Vector4(x1_x0.x, x1_x0.y, x1_x0.z, 0));
        shapeMatrix.SetRow(1, new Vector4(x2_x0.x, x2_x0.y, x2_x0.z, 0));
        shapeMatrix.SetRow(2, new Vector4(x3_x0.x, x3_x0.y, x3_x0.z, 0));

        return shapeMatrix;
    }

    public float CalculateC(Triangle triangle)
    {
        Matrix4x4 currentMatrix = CalculateShapeMatrix(triangle);

        // 计算两个3x3矩阵的行列式差值
        float detCurrent = currentMatrix.m00 * (currentMatrix.m11 * currentMatrix.m22 - currentMatrix.m12 * currentMatrix.m21) -
                           currentMatrix.m01 * (currentMatrix.m10 * currentMatrix.m22 - currentMatrix.m12 * currentMatrix.m20) +
                           currentMatrix.m02 * (currentMatrix.m10 * currentMatrix.m21 - currentMatrix.m11 * currentMatrix.m20);

        float detInitial = triangle.initialMatrix.m00 * (triangle.initialMatrix.m11 * triangle.initialMatrix.m22 - triangle.initialMatrix.m12 * triangle.initialMatrix.m21) -
                           triangle.initialMatrix.m01 * (triangle.initialMatrix.m10 * triangle.initialMatrix.m22 - triangle.initialMatrix.m12 * triangle.initialMatrix.m20) +
                           triangle.initialMatrix.m02 * (triangle.initialMatrix.m10 * triangle.initialMatrix.m21 - triangle.initialMatrix.m11 * triangle.initialMatrix.m20);

        return CalcDet(currentMatrix) - CalcDet(triangle.initialMatrix);//detCurrent - detInitial;
    }

    public float CalcDet(Matrix4x4 mat)
    {
        Vector3 a = mat.GetRow(0);
        Vector3 b = mat.GetRow(1);
        Vector3 c = mat.GetRow(2);

        return Vector3.Dot(a, Vector3.Cross(b, c));
    }

    public void CalculateGrad(Triangle triangle, out Vector3 grad1, out Vector3 grad2, out Vector3 grad3, out Vector3 grad0)
    {
        Vector3 curCentroid = CalculateCentroid(triangle);
        Vector3 p1 = mySolver.pos[triangle.pointsId[0]] * 100;
        Vector3 p2 = mySolver.pos[triangle.pointsId[1]] * 100;
        Vector3 p3 = mySolver.pos[triangle.pointsId[2]] * 100;

        // 计算各个点相对质心的向量
        Vector3 p1_p0 = p1 - curCentroid;  // x1 - x0
        Vector3 p2_p0 = p2 - curCentroid;  // x2 - x0
        Vector3 p3_p0 = p3 - curCentroid;  // x3 - x0

        // 根据公式计算各个点的梯度
        grad1 = Vector3.Cross(p2_p0, p3_p0);  // ∇x1C = (x2 - x0) × (x3 - x0)
        grad2 = Vector3.Cross(p3_p0, p1_p0);  // ∇x2C = (x3 - x0) × (x1 - x0)
        grad3 = Vector3.Cross(p1_p0, p2_p0);  // ∇x3C = (x1 - x0) × (x2 - x0)

        // 计算质心的梯度
        grad0 = -(grad1 + grad2 + grad3);     // ∇x0C = -∇x1C - ∇x2C - ∇x3C
    }

    public override void SolveConstraint(float dt)
    {
        float compliance = stiff / (dt * dt);

        foreach (var triangle in mySolver.triangles)
        {
            // 获取顶点位置

            // 计算约束值
            float C = CalculateC(triangle);
            //Debug.Log("约束至"+C);

            // 计算梯度
            Vector3 grad1, grad2, grad3, grad0;
            CalculateGrad(triangle, out grad1, out grad2, out grad3, out grad0);

            // 计算分母
            float w1 = 1; //mySolver.invMass[triangle.pointsId[0]];
            float w2 = 1;// mySolver.invMass[triangle.pointsId[1]];
            float w3 = 1;// mySolver.invMass[triangle.pointsId[2]];

            float denominator =
                w1 * grad1.magnitude * grad1.magnitude +
                w2 * grad2.magnitude * grad2.magnitude +
                w3 * grad3.magnitude * grad3.magnitude +
                compliance;

            // 计算delta_lambda并累加到总lambda
            float deltaLambda = denominator != 0.0f ? -(C + compliance * triangle.lambda) / denominator : 0.0f;
            triangle.lambda += deltaLambda;
            //Debug.Log("lambda:" + triangle.lambda);

            //Debug.Log($"{w1} * {deltaLambda} * {grad1}");
            // 使用累加的lambda更新位置
            if (w1 > 0) mySolver.pos[triangle.pointsId[0]] += w1 * deltaLambda * grad1 / 100;
            if (w2 > 0) mySolver.pos[triangle.pointsId[1]] += w2 * deltaLambda * grad2 / 100;
            if (w3 > 0) mySolver.pos[triangle.pointsId[2]] += w3 * deltaLambda * grad3 / 100;
        }
    }
}

public class STVKElasticConstraint : Constraint
{
    Dictionary<KeyValuePair<int, int>, float> distances = new Dictionary<KeyValuePair<int, int>, float>();

    int[] pointAs;
    int[] pointBs;
    float[] distance;

    float[] lambdas;

    float alpha;
    public STVKElasticConstraint(ClothXPBDSolver solver) : base(solver)
    {
        for (int i = 0; i < mySolver.triangles.Count(); i++)
        {
            List<int> points = mySolver.triangles[i].pointsId.ToList();
            points.Sort();
            var pair = points.SelectMany((first, index) => points.Skip(index + 1), (first, second) => new { first, second });
            foreach (var point in pair)
            {
                if (!distances.ContainsKey(new(point.first, point.second)))
                {
                    distances.Add(new(point.first, point.second), (mySolver.pos[point.first] - mySolver.pos[point.second]).magnitude);
                }
            }
        }

        pointAs = new int[distances.Count];
        pointBs = new int[distances.Count];
        distance = new float[distances.Count];
        int j = 0;
        foreach (var pair in distances)
        {
            pointAs[j] = pair.Key.Key;
            pointBs[j] = pair.Key.Value;
            distance[j] = pair.Value;
            j++;
        }

        lambdas = new float[distances.Count];
    }

    public override void SolveConstraint(float dt)
    {
        stiff = 1;
        alpha = stiff / (math.pow(dt, 2));

        for (int i = 0; i < distance.Length; i++)
        {
            int particleA = pointAs[i];
            int particleB = pointBs[i];

            // 计算约束梯度
            Vector3 gradientC = CalculateGradC(particleA, particleB, i);
            float gradientMagnitude = gradientC.magnitude == 0 ? 0.1f : gradientC.magnitude;

            // 计算拉格朗日乘子的增量
            float deltaLambda = CalcDeltaLambda(particleA, particleB, gradientMagnitude, i);

            // 计算位置修正
            Vector3 positionCorrection = gradientC * deltaLambda;

            // 应用位置修正
            mySolver.pos[particleB] += positionCorrection * mySolver.invMass[particleB];
            mySolver.pos[particleA] -= positionCorrection * mySolver.invMass[particleA];

            lambdas[i] += deltaLambda;
        }
    }

    public float CalculateC(int particleA, int particleB, int i)
    {
        Vector3 displacement = mySolver.pos[particleA] - mySolver.pos[particleB];
        float currentDistance = displacement.magnitude;
        float restDistance = distance[i];

        float distanceError = math.pow(currentDistance, 2) - math.pow(restDistance, 2);
        float exponentialTerm = math.exp(0.5f * math.pow(distanceError, 2));

        return stiff * math.pow(1 - exponentialTerm, 2);
    }

    public Vector3 CalculateGradC(int particleA, int particleB, int i)
    {
        Vector3 displacement = mySolver.pos[particleA] - mySolver.pos[particleB];
        float restDistance = distance[i];

        // 计算距离误差
        float distanceError = math.pow(displacement.magnitude, 2) - math.pow(restDistance, 2);
        float exponentialTerm = math.exp(0.5f * math.pow(distanceError, 2));

        // 计算梯度各个部分
        float stiffnessTerm = 2 * stiff * (1 - exponentialTerm);
        float exponentialFactor = -exponentialTerm;
        Vector3 directionTerm = 2 * distanceError * (mySolver.pos[particleB] - mySolver.pos[particleA]);

        return stiffnessTerm * exponentialFactor * directionTerm;
    }

    public float CalcDeltaLambda(int particleA, int particleB, float gradientMagnitude, int i)
    {
        float constraintValue = CalculateC(particleA, particleB, i);
        float dampingTerm = alpha * lambdas[i];
        float denominator = gradientMagnitude * gradientMagnitude *
                           (mySolver.invMass[particleA] + mySolver.invMass[particleB]) + alpha;

        return -(constraintValue + dampingTerm) / denominator;
    }

    public override void ResetLambda()
    {
        lambdas = new float[distances.Count];
    }
}

public class ElasticConstraint : Constraint
{

    int[] pointAs;
    int[] pointBs;
    float[] distance;

    float[] lambdas;

    float alpha;
    public ElasticConstraint(ClothXPBDSolver solver) : base(solver)
    {
        stiff = 0f;
        Dictionary<KeyValuePair<int, int>, float> distances = new Dictionary<KeyValuePair<int, int>, float>();

        for (int i = 0; i < mySolver.triangles.Count(); i++)
        {
            List<int> points = mySolver.triangles[i].pointsId.ToList();
            points.Sort();
            var pair = points.SelectMany((first, index) => points.Skip(index + 1), (first, second) => new { first, second });
            foreach (var point in pair)
            {
                if (!distances.ContainsKey(new(point.first, point.second)))
                {
                    distances.Add(new(point.first, point.second), (mySolver.pos[point.first] - mySolver.pos[point.second]).magnitude);
                }
            }
        }

        pointAs = new int[distances.Count];
        pointBs = new int[distances.Count];
        distance = new float[distances.Count];
        int j = 0;
        foreach (var pair in distances)
        {
            pointAs[j] = pair.Key.Key;
            pointBs[j] = pair.Key.Value;
            distance[j] = pair.Value;
            j++;
        }

        lambdas = new float[distances.Count];
    }

    public override void SolveConstraint(float dt)
    {
        alpha = stiff / (math.pow(dt, 2));

        for (int i = 0; i < distance.Length; i++)
        {
            int particleA = pointAs[i];
            int particleB = pointBs[i];

            // 计算约束梯度
            Vector3 gradientC = CalculateGradC(particleA, particleB, i);
            float gradientMagnitude = gradientC.magnitude == 0 ? 0.1f : gradientC.magnitude;

            // 计算拉格朗日乘子的增量
            float deltaLambda = CalcDeltaLambda(particleA, particleB, gradientMagnitude, i);

            // 计算位置修正
            Vector3 positionCorrection = gradientC * deltaLambda;

            //Debug.Log(positionCorrection * mySolver.invMass[particleB]);
            // 应用位置修正
            mySolver.pos[particleB] -= positionCorrection * mySolver.invMass[particleB];
            mySolver.pos[particleA] += positionCorrection * mySolver.invMass[particleA];

            lambdas[i] += deltaLambda;
        }
    }

    public float CalculateC(int particleA, int particleB, int i)
    {
        Vector3 displacement = mySolver.pos[particleA] - mySolver.pos[particleB];
        float currentDistance = displacement.magnitude;
        float restDistance = distance[i];

        float distanceError = currentDistance - restDistance;
        float exponentialTerm = math.exp(distanceError);

        return stiff * math.pow(1 - exponentialTerm, 2);
    }

    public Vector3 CalculateGradC(int particleA, int particleB, int i)
    {
        Vector3 displacement = mySolver.pos[particleA] - mySolver.pos[particleB];
        float restDistance = distance[i];

        // 计算距离误差
        float distanceError = displacement.magnitude - restDistance;
        float exponentialTerm = math.exp(distanceError);

        // 计算梯度各个部分
        float stiffnessTerm = 2 * stiff * (1 - exponentialTerm);
        float exponentialFactor = -exponentialTerm;
        Vector3 directionTerm = (1 / displacement.magnitude) * displacement;

        return stiffnessTerm * exponentialFactor * directionTerm;
    }

    public float CalcDeltaLambda(int particleA, int particleB, float gradientMagnitude, int i)
    {
        float constraintValue = CalculateC(particleA, particleB, i);
        float dampingTerm = alpha * lambdas[i];
        float denominator = gradientMagnitude * gradientMagnitude *
                           (mySolver.invMass[particleA] + mySolver.invMass[particleB]) + alpha;

        return -(constraintValue + dampingTerm) / denominator;
    }

    public override void ResetLambda()
    {
        lambdas = new float[distance.Length];
    }
}

public class Volume2Constraint : Constraint
{
    float lambda = 0f;

    float stiff = 0.9f;

    float oriVolume = 0f;
    Vector3[] curNormals;

    public Volume2Constraint(XPBDSolver solver) : base(solver)
    {
        // 预先计算每个顶点的“面积加权法线”（只需一次）
        // 以便后续计算体积和梯度。
        curNormals = new Vector3[mySolver.numParticles];
        oriVolume = 0;
        for (int t = 0; t < mySolver.triangles.Length; t++)
        {
            Triangle tri = mySolver.triangles[t];
            Vector3 p0 = mySolver.pos[tri.pointsId[0]];
            Vector3 p1 = mySolver.pos[tri.pointsId[1]];
            Vector3 p2 = mySolver.pos[tri.pointsId[2]];

            // 计算面积和法线
            Vector3 faceNormal = Vector3.Cross(p1 - p0, p2 - p0);
            float area = faceNormal.magnitude * 0.5f;
            Vector3 normal = faceNormal.normalized;

            // 计算三个顶点位置的和与法线的点积
            Vector3 posSum = p0 + p1 + p2;
            oriVolume += (1.0f / 3.0f) * area * Vector3.Dot(posSum, normal);
        }
        //Debug.Log("ori"+oriVolume);
    }

    public override void ResetLambda()
    {
        lambda = 0f;
    }

    void ComputeAreaWeightedNormals(Vector3[] v)
    {
        // 初始化
        for (int i = 0; i < mySolver.numParticles; i++)
        {
            v[i] = Vector3.zero;
        }

        // 对于每个三角形，把三角形对应的 “面积*法线” 加到它的三个顶点上
        for (int t = 0; t < mySolver.triangles.Length; t++)
        {
            Triangle tri = mySolver.triangles[t];
            int i0 = tri.pointsId[0];
            int i1 = tri.pointsId[1];
            int i2 = tri.pointsId[2];

            // 初始顶点位置
            Vector3 p0 = mySolver.pos[i0];
            Vector3 p1 = mySolver.pos[i1];
            Vector3 p2 = mySolver.pos[i2];

            // 三角形法线 (p1 - p0) x (p2 - p0)
            Vector3 faceNormal = Vector3.Cross(p1 - p0, p2 - p0);

            // 面积的一半
            float area = faceNormal.magnitude * 0.5f;
            // 归一化法线
            Vector3 normal = faceNormal.normalized;

            // 面积加权法线
            Vector3 areaWeighted = normal * area;

            v[i0] += areaWeighted;
            v[i1] += areaWeighted;
            v[i2] += areaWeighted;
        }
    }

    public override void SolveConstraint(float dt)
    {
        ComputeAreaWeightedNormals(curNormals);
        // 1) 计算当前体积
        float currentVolume = 0;
        for (int t = 0; t < mySolver.triangles.Length; t++)
        {
            Triangle tri = mySolver.triangles[t];
            Vector3 p0 = mySolver.pos[tri.pointsId[0]];
            Vector3 p1 = mySolver.pos[tri.pointsId[1]];
            Vector3 p2 = mySolver.pos[tri.pointsId[2]];

            // 计算面积和法线
            Vector3 faceNormal = Vector3.Cross(p1 - p0, p2 - p0);
            float area = faceNormal.magnitude * 0.5f;
            Vector3 normal = faceNormal.normalized;

            // 计算三个顶点位置的和与法线的点积
            Vector3 posSum = p0 + p1 + p2;
            currentVolume += (1.0f / 3.0f) * area * Vector3.Dot(posSum, normal);
        }
        //Debug.Log("cur" + currentVolume);

        float C = currentVolume - oriVolume;

        // 如果你想要只在物体“被压缩”时才修正，也可以加个 if(C < 0) 才纠正等逻辑
        // 但一般全程都修正，以保证体积不发散。

        // 3) 计算 XPBD 中的 α = stiffness / dt^2
        float alpha = stiff / (dt * dt);


        // 先把 gradC 和 sumWGradCSq 计算好
        Vector3[] gradC = new Vector3[mySolver.numParticles];
        float sumWGradCSq = 0f;
        for (int i = 0; i < mySolver.numParticles; i++)
        {
            // 质量为 1/invMass[i]
            float w_i = mySolver.invMass[i];  // 注意：这里是“质量的倒数”

            if (Mathf.Abs(w_i) < 1e-10f)
            {
                gradC[i] = Vector3.zero;
                continue; // 跳过不可动点
            }

            gradC[i] = curNormals[i] / 3f;

            // |gradC[i]|^2
            float gradSq = gradC[i].sqrMagnitude;

            sumWGradCSq += w_i * gradSq;
        }

        // 5) XPBD 的拉格朗日乘子增量
        //   deltalambda = -(C + alpha * lambda) / ( sumWGradCSq + alpha )
        float deltalambda = -(C + alpha * lambda) / (sumWGradCSq + alpha);

        // 更新本约束的 λ
        lambda += deltalambda;

        // 6) 根据 deltalambda 给每个顶点做位置修正
        //   Δxᵢ = - (deltalambda * wᵢ) * gradC[i]
        //   (注意符号，有些实现里写成 “+”，主要看公式里 C 的正负，这里习惯用“负号”)
        for (int i = 0; i < mySolver.numParticles; i++)
        {
            float w_i = mySolver.invMass[i];
            if (Mathf.Abs(w_i) < 1e-10f)
                continue;

            // XPBD 位置修正
            Vector3 deltaX = deltalambda * w_i * gradC[i];
            mySolver.pos[i] += deltaX;
        }

        // 至此，体积约束的迭代就完成一次。
    }
}

public class BendingConstraint : Constraint
{
    int[] x1;
    int[] x2;
    int[] x3;
    int[] x4;

    Vector3[] initialN1;
    Vector3[] initialN2;

    float[] lambdas;
    float alpha;
    public BendingConstraint(ClothXPBDSolver solver) : base(solver)
    {
        stiff = 0;
        Dictionary<KeyValuePair<int, int>, List<int>> triangleEdge = new Dictionary<KeyValuePair<int, int>, List<int>>();

        for (int i = 0; i < mySolver.triangles.Count(); i++)
        {
            List<int> points = mySolver.triangles[i].pointsId.ToList();
            points.Sort();
            var pair = points.SelectMany((first, index) => points.Skip(index + 1), (first, second) => new { first, second });
            foreach (var point in pair)
            {
                if (!triangleEdge.ContainsKey(new(point.first, point.second)))
                {
                    triangleEdge.Add(new(point.first, point.second), new List<int>() { i });
                }
                else
                {
                    triangleEdge[new(point.first, point.second)].Add(i);
                }
            }
        }

        List<int> x1 = new List<int>();
        List<int> x2 = new List<int>();
        List<int> x3 = new List<int>();
        List<int> x4 = new List<int>();
        List<Vector3> n1 = new List<Vector3>();
        List<Vector3> n2 = new List<Vector3>();

        foreach (var item in triangleEdge)
        {
            var pairTriangles = item.Value.SelectMany((first, index) => item.Value.Skip(index + 1), (first, second) => new { first, second });
            foreach (var triangle in pairTriangles)
            {
                x1.Add(item.Key.Key);
                x2.Add(item.Key.Value);
                x3.Add(mySolver.triangles[triangle.first].pointsId.Where(p => p != x1.Last() && p != x2.Last()).First());
                x4.Add(mySolver.triangles[triangle.second].pointsId.Where(p => p != x1.Last() && p != x2.Last()).First());

                Vector3 p2 = mySolver.pos[x2.Last()] - mySolver.pos[x1.Last()];
                Vector3 p3 = mySolver.pos[x3.Last()] - mySolver.pos[x1.Last()];
                Vector3 p4 = mySolver.pos[x4.Last()] - mySolver.pos[x1.Last()];
                n1.Add(CalcNormal(Vector3.zero, p2, p3));
                n2.Add(CalcNormal(Vector3.zero, p2, p4));

                if (1 - math.abs(Vector3.Dot(n1.Last(), n2.Last())) <= 0.001f)
                {
                    x1.RemoveAt(x1.Count - 1);
                    x2.RemoveAt(x1.Count);
                    x3.RemoveAt(x1.Count);
                    x4.RemoveAt(x1.Count);
                    n1.RemoveAt(x1.Count);
                    n2.RemoveAt(x1.Count);
                }

                //Debug.Log($"{x1.Last()} + {x2.Last()} + {mySolver.triangles[triangle.first].pointsId[0]} + {mySolver.triangles[triangle.first].pointsId[1]} + {mySolver.triangles[triangle.first].pointsId[2]} + {mySolver.triangles[triangle.second].pointsId[0]} + {mySolver.triangles[triangle.second].pointsId[1]} + {mySolver.triangles[triangle.second].pointsId[2]}");

            }
        }

        this.x1 = x1.ToArray();
        this.x2 = x2.ToArray();
        this.x3 = x3.ToArray();
        this.x4 = x4.ToArray();
        initialN1 = n1.ToArray();
        initialN2 = n2.ToArray();

        lambdas = new float[x1.Count];
    }

    public Vector3 CalcNormal(Vector3 x1, Vector3 x2, Vector3 x3)
    {
        Vector3 normalized = Vector3.Cross(x2 - x1, x3 - x1).normalized;

        return CheckValid(normalized);
    }

    public Vector3 CheckValid(Vector3 vec)
    {
        if (float.IsNaN(vec.x) || float.IsNaN(vec.y) || float.IsNaN(vec.z) ||
        float.IsInfinity(vec.x) || float.IsInfinity(vec.y) || float.IsInfinity(vec.z))
        {
            Debug.Log($"wuwawa!!!");
            return Vector3.zero;
        }
        return vec;
    }

    public override void SolveConstraint(float dt)
    {
        alpha = stiff / (math.pow(dt, 2));

        for (int i = 0; i < x1.Length; i++)
        {
            Vector3 p2 = mySolver.pos[x2[i]] - mySolver.pos[x1[i]];
            Vector3 p3 = mySolver.pos[x3[i]] - mySolver.pos[x1[i]];
            Vector3 p4 = mySolver.pos[x4[i]] - mySolver.pos[x1[i]];

            Vector3 n1 = CalcNormal(Vector3.zero, p2, p3);
            Vector3 n2 = CalcNormal(Vector3.zero, p2, p4);
            float d = Vector3.Dot(n1, n2);
            //Debug.Log($"{mySolver.pos[x1[i]]} {p2} {p3} {p4} {n1} {n2} {d}");
            float C = CalculateC(n1, n2, d, i);

            (Vector3 grad1, Vector3 grad2, Vector3 grad3, Vector3 grad4) = CalculateGradC(p2, p3, p4, n1, n2, d, i);

            grad1 = CheckValid(grad1);
            grad2 = CheckValid(grad2);
            grad3 = CheckValid(grad3);
            grad4 = CheckValid(grad4);

            //Debug.Log($"{grad1} {grad2} {grad3} {grad4} {d}");
            float deltaLambda = -(C * math.sqrt(1 - d * d) + alpha * lambdas[i]) / (mySolver.invMass[x1[i]] * grad1.magnitude * grad1.magnitude + mySolver.invMass[x2[i]] * grad2.magnitude * grad2.magnitude + mySolver.invMass[x3[i]] * grad3.magnitude * grad3.magnitude + mySolver.invMass[x4[i]] * grad4.magnitude * grad4.magnitude + alpha);

            lambdas[i] += deltaLambda;
            //Debug.Log(mySolver.invMass[x1[i]] * grad1 * deltaLambda);

            mySolver.pos[x1[i]] += mySolver.invMass[x1[i]] * grad1 * deltaLambda;
            mySolver.pos[x2[i]] += mySolver.invMass[x2[i]] * grad2 * deltaLambda;
            mySolver.pos[x3[i]] += mySolver.invMass[x3[i]] * grad3 * deltaLambda;
            mySolver.pos[x4[i]] += mySolver.invMass[x4[i]] * grad4 * deltaLambda;
        }
    }

    public float CalculateC(Vector3 n1, Vector3 n2, float d, int i)
    {
        return math.acos(d) - math.acos(Vector3.Dot(initialN1[i], initialN2[i]));
    }

    public (Vector3, Vector3, Vector3, Vector3) CalculateGradC(Vector3 p2, Vector3 p3, Vector3 p4, Vector3 n1, Vector3 n2, float d, int i)
    {
        Vector3 q3 = Vector3.Cross(p2, p3).magnitude != 0 ? (Vector3.Cross(p2, n2) + (Vector3.Cross(n1, p2) * d)) / Vector3.Cross(p2, p3).magnitude : Vector3.zero;

        Vector3 q4 = Vector3.Cross(p2, p4).magnitude != 0 ? (Vector3.Cross(p2, n1) + (Vector3.Cross(n2, p2) * d)) / Vector3.Cross(p2, p4).magnitude : Vector3.zero;

        Vector3 q2 = (Vector3.Cross(p2, p3).magnitude != 0 && Vector3.Cross(p2, p4).magnitude != 0) ? -(Vector3.Cross(p3, n2) + (Vector3.Cross(n1, p3) * d)) / Vector3.Cross(p2, p3).magnitude - (Vector3.Cross(p4, n1) + (Vector3.Cross(n2, p4) * d)) / Vector3.Cross(p2, p4).magnitude : Vector3.zero;

        if (math.abs(d - 1) < 0.001f) d = 0.9f;
        float down = 1;/// math.sqrt(1 - d * d);
        //Debug.Log($"{q2} {q3} {q4}");
        return ((-q2 - q3 - q4) * down, q2 * down, q3 * down, q4 * down);
    }

    public override void ResetLambda()
    {
        lambdas = new float[x1.Count()];
    }
}

public class SectionVolume2Constraint : Constraint
{
    // 存储所有section对的信息
    private class SectionPair
    {
        public int[] section1;
        public int[] section2;
        public HashSet<int> vertices;
        public float lambda;
        public float oriVolume;
        public Vector3 ghostPoint1; // section1的圆心幽灵点
        public Vector3 ghostPoint2; // section2的圆心幽灵点
        public Vector3 sectionDirection;
    }

    private List<SectionPair> sectionPairs;
    private float stiff = 0.5f;
    private Vector3[] curNormals;


    private Vector3 CalculateSectionCenter(int[] sectionVertices)
    {
        Vector3 center = Vector3.zero;
        foreach (int idx in sectionVertices)
        {
            center += mySolver.pos[idx];
        }
        return center / sectionVertices.Length;
    }

    public SectionVolume2Constraint(ClothXPBDSolver solver) : base(solver)
    {
        sectionPairs = new List<SectionPair>();
        curNormals = new Vector3[mySolver.numParticles];

        // 创建所有section对
        if (((ClothXPBDSolver)mySolver).sections.Length > 0)
        {
            // 第一段：端点0和第一个section
            CreateSectionPair(new int[] { 0 }, ((ClothXPBDSolver)mySolver).sections[0]);

            // 中间sections之间
            for (int i = 0; i < ((ClothXPBDSolver)mySolver).sections.Length - 1; i++)
            {
                CreateSectionPair(((ClothXPBDSolver)mySolver).sections[i], ((ClothXPBDSolver)mySolver).sections[i + 1]);
            }

            // 最后一段：最后一个section和末端点
            CreateSectionPair(((ClothXPBDSolver)mySolver).sections[((ClothXPBDSolver)mySolver).sections.Length - 1],
                new int[] { mySolver.numParticles - 1 });
        }
    }

    public SectionVolume2Constraint(XPBDSolver solver) : base(solver)
    {
    }

    private void CreateSectionPair(int[] firstSection, int[] secondSection)
    {
        SectionPair pair = new SectionPair
        {
            section1 = firstSection,
            section2 = secondSection,
            vertices = new HashSet<int>(),
            lambda = 0f
        };

        // 合并顶点
        foreach (int v in firstSection) pair.vertices.Add(v);
        foreach (int v in secondSection) pair.vertices.Add(v);

        // 计算两个section的幽灵点（圆心）
        pair.ghostPoint1 = CalculateSectionCenter(firstSection);
        pair.ghostPoint2 = CalculateSectionCenter(secondSection);

        // 计算section的方向（从section1指向section2）
        pair.sectionDirection = (pair.ghostPoint2 - pair.ghostPoint1).normalized;

        // 确保section点的顺序一致性
        if (firstSection.Length > 1)
        {
            EnsureSectionOrientation(firstSection, pair.ghostPoint1, pair.sectionDirection);
        }
        if (secondSection.Length > 1)
        {
            EnsureSectionOrientation(secondSection, pair.ghostPoint2, pair.sectionDirection);
        }

        pair.oriVolume = CalculateVolumeWithGhostPoints(pair);
        sectionPairs.Add(pair);
    }

    private void EnsureSectionOrientation(int[] sectionPoints, Vector3 ghostPoint, Vector3 sectionDir)
    {
        if (sectionPoints.Length <= 2) return;

        // 计算section的法线方向
        Vector3 firstEdge = mySolver.pos[sectionPoints[1]] - mySolver.pos[sectionPoints[0]];
        Vector3 center = ghostPoint;
        Vector3 toCenter = (center - mySolver.pos[sectionPoints[0]]).normalized;
        Vector3 currentNormal = Vector3.Cross(firstEdge, toCenter);

        // 如果法线方向不一致，翻转点的顺序
        if (Vector3.Dot(currentNormal, sectionDir) < 0)
        {
            System.Array.Reverse(sectionPoints);
        }
    }


    private float CalculateVolumeWithGhostPoints(SectionPair pair)
    {
        float volume = 0f;

        // 1. 计算原有三角形的体积贡献
        for (int t = 0; t < mySolver.triangles.Length; t++)
        {
            Triangle tri = mySolver.triangles[t];
            if (!IsTriangleInSection(tri, pair.vertices)) continue;

            volume += CalculateTriangleVolume(
                mySolver.pos[tri.pointsId[0]],
                mySolver.pos[tri.pointsId[1]],
                mySolver.pos[tri.pointsId[2]],
                false,
                pair
            );
        }

        // 2. 计算section1的封闭面（确保逆时针顺序，法线朝外）
        if (pair.section1.Length > 1)
        {
            for (int i = 0; i < pair.section1.Length; i++)
            {
                int current = pair.section1[i];
                int next = pair.section1[(i + 1) % pair.section1.Length];

                volume += CalculateTriangleVolume(
                    mySolver.pos[current],
                    mySolver.pos[next],
                    pair.ghostPoint1,
                    true,
                    pair
                );
            }
        }

        // 3. 计算section2的封闭面（确保顺时针顺序，法线朝内）
        if (pair.section2.Length > 1)
        {
            for (int i = 0; i < pair.section2.Length; i++)
            {
                int current = pair.section2[i];
                int next = pair.section2[(i + 1) % pair.section2.Length];

                volume += CalculateTriangleVolume(
                    mySolver.pos[current],
                    mySolver.pos[next],
                    pair.ghostPoint2,
                    false,
                    pair
                );
            }
        }

        return volume;
    }

    private float CalculateTriangleVolume(Vector3 p0, Vector3 p1, Vector3 p2, bool flipNormal, SectionPair pair)
    {
        // 确保点的顺序产生正确的法线方向
        Vector3 edge1 = p1 - p0;
        Vector3 edge2 = p2 - p0;
        Vector3 faceNormal = Vector3.Cross(edge1, edge2);

        if (flipNormal)
        {
            faceNormal = -faceNormal;
        }

        float area = faceNormal.magnitude * 0.5f;
        Vector3 normal = faceNormal.normalized;

        // 不再需要额外的法线方向检查，因为我们在创建时已经确保了正确的方向
        Vector3 posSum = p0 + p1 + p2;
        return (1.0f / 9.0f) * area * Vector3.Dot(posSum, normal);
    }



    private bool IsTriangleInSection(Triangle tri, HashSet<int> sectionVertices)
    {
        int verticesInSection = 0;
        if (sectionVertices.Contains(tri.pointsId[0])) verticesInSection++;
        if (sectionVertices.Contains(tri.pointsId[1])) verticesInSection++;
        if (sectionVertices.Contains(tri.pointsId[2])) verticesInSection++;
        return verticesInSection >= 2;
    }

    public override void ResetLambda()
    {
        foreach (var pair in sectionPairs)
        {
            pair.lambda = 0f;
        }
    }

    void ComputeAreaWeightedNormals(Vector3[] v)
    {
        // 初始化所有法线为零
        for (int i = 0; i < mySolver.numParticles; i++)
        {
            v[i] = Vector3.zero;
        }

        // 对每个section对计算法线
        foreach (var pair in sectionPairs)
        {
            for (int t = 0; t < mySolver.triangles.Length; t++)
            {
                Triangle tri = mySolver.triangles[t];
                if (!IsTriangleInSection(tri, pair.vertices)) continue;

                int i0 = tri.pointsId[0];
                int i1 = tri.pointsId[1];
                int i2 = tri.pointsId[2];

                Vector3 p0 = mySolver.pos[i0];
                Vector3 p1 = mySolver.pos[i1];
                Vector3 p2 = mySolver.pos[i2];

                Vector3 faceNormal = Vector3.Cross(p1 - p0, p2 - p0);
                float area = faceNormal.magnitude * 0.5f;
                Vector3 normal = faceNormal.normalized;
                Vector3 areaWeighted = normal * area;

                v[i0] += areaWeighted;
                v[i1] += areaWeighted;
                v[i2] += areaWeighted;
            }
        }
    }

    public override void SolveConstraint(float dt)
    {
        ComputeAreaWeightedNormals(curNormals);

        foreach (var pair in sectionPairs)
        {
            // 更新幽灵点位置
            pair.ghostPoint1 = CalculateSectionCenter(pair.section1);
            pair.ghostPoint2 = CalculateSectionCenter(pair.section2);

            float currentVolume = CalculateVolumeWithGhostPoints(pair);
            float C = currentVolume - pair.oriVolume;
            float alpha = stiff / (dt * dt);

            Vector3[] gradC = new Vector3[mySolver.numParticles];
            Vector3 gradC_ghost1 = Vector3.zero;  // 幽灵点1的梯度
            Vector3 gradC_ghost2 = Vector3.zero;  // 幽灵点2的梯度
            float sumWGradCSq = 0f;

            // 计算实际顶点的梯度
            foreach (int i in pair.vertices)
            {
                gradC[i] = curNormals[i] / 3f;

                float w_i = mySolver.invMass[i];
                if (w_i > 1e-10f)
                {
                    sumWGradCSq += w_i * gradC[i].sqrMagnitude;
                }
            }

            // 计算幽灵点的梯度贡献
            // Section1的封闭面对幽灵点1的梯度
            for (int i = 0; i < pair.section1.Length; i++)
            {
                int current = pair.section1[i];
                int next = pair.section1[(i + 1) % pair.section1.Length];
                Vector3 p0 = mySolver.pos[current];
                Vector3 p1 = mySolver.pos[next];

                Vector3 edge = p1 - p0;
                Vector3 normal = Vector3.Cross(edge, (p0 - pair.ghostPoint1).normalized);
                gradC_ghost1 += normal / 3f;
            }

            // Section2的封闭面对幽灵点2的梯度
            for (int i = 0; i < pair.section2.Length; i++)
            {
                int current = pair.section2[i];
                int next = pair.section2[(i + 1) % pair.section2.Length];
                Vector3 p0 = mySolver.pos[current];
                Vector3 p1 = mySolver.pos[next];

                Vector3 edge = p1 - p0;
                Vector3 normal = Vector3.Cross(edge, (p0 - pair.ghostPoint2).normalized);
                gradC_ghost2 += normal / 3f;
            }

            // 假设幽灵点的质量为实际顶点平均质量
            float avgMass = 0f;
            int validPoints = 0;
            foreach (int i in pair.vertices)
            {
                if (mySolver.invMass[i] > 1e-10f)
                {
                    avgMass += 1f / mySolver.invMass[i];
                    validPoints++;
                }
            }
            if (validPoints > 0)
            {
                avgMass /= validPoints;
                float ghostInvMass = 1f / avgMass;

                // 将幽灵点的贡献添加到sumWGradCSq
                sumWGradCSq += ghostInvMass * gradC_ghost1.sqrMagnitude;
                sumWGradCSq += ghostInvMass * gradC_ghost2.sqrMagnitude;
            }

            // 计算位置修正
            float deltalambda = -(C + alpha * pair.lambda) / (sumWGradCSq + alpha);
            pair.lambda += deltalambda;

            // 更新实际顶点的位置
            foreach (int i in pair.vertices)
            {
                float w_i = mySolver.invMass[i];
                if (w_i < 1e-10f)
                    continue;

                Vector3 deltaX = deltalambda * w_i * gradC[i];
                mySolver.pos[i] += deltaX;
            }

            // 幽灵点不需要更新位置，但它们的梯度已经影响了约束求解
        }
    }
}

public class SectionAreaConstriant : Constraint
{
    int[] x1;
    int[] x2;
    int[] x3;
    Vector3[] ori2_1;
    Vector3[] ori3_1;

    float alpha;
    float[] lambdas;
    public SectionAreaConstriant(ClothXPBDSolver solver) : base(solver)
    {
        List<int> x1 = new List<int>();
        List<int> x2 = new List<int>();
        List<int> x3 = new List<int>();
        List<Vector3> x2_1 = new List<Vector3>();
        List<Vector3> x3_1 = new List<Vector3>();

        for (int i = 0; i < solver.sections.Count(); i++)
        {
            for (int j = 0; j < solver.subdivision - 2; j++)
            {
                x1.Add(solver.sections[i][0]);
                x2.Add(solver.sections[i][j + 1]);
                x3.Add(solver.sections[i][j + 2]);
                x2_1.Add(solver.pos[x2.Last()] - solver.pos[x1.Last()]);
                x3_1.Add(solver.pos[x3.Last()] - solver.pos[x1.Last()]);
                //Debug.Log($"{x1.Last()}, {x2.Last()}, {x3.Last()}");
            }
        }

        this.x1 = x1.ToArray();
        this.x2 = x2.ToArray();
        this.x3 = x3.ToArray();
        this.ori2_1 = x2_1.ToArray();
        this.ori3_1 = x3_1.ToArray();

        lambdas = new float[x1.Count()];
    }

    public override void ResetLambda()
    {
        lambdas = new float[x1.Count()];
    }

    public override void SolveConstraint(float dt)
    {
        alpha = stiff / (dt * dt);

        for (int i = 0; i < x1.Count(); i++)
        {
            Vector3 x2_x1 = mySolver.pos[x2[i]] * 20 - mySolver.pos[x1[i]] * 20;
            Vector3 x3_x1 = mySolver.pos[x3[i]] * 20 - mySolver.pos[x1[i]] * 20;
            // Debug.Log($"{x2_x1}, {x3_x1}");
            float areaConstraint = math.pow(Vector3.Cross(x2_x1, x3_x1).magnitude, 2) - math.pow(Vector3.Cross(ori2_1[i] * 20, ori3_1[i] * 20).magnitude, 2);
            //Debug.Log($"{areaConstraint}, {Vector3.Cross(x2_x1, x3_x1)}");
            Vector3 grad2 = Vector3.Cross(x3_x1, Vector3.Cross(x2_x1, x3_x1)) * 2;
            Vector3 grad3 = Vector3.Cross(x2_x1, Vector3.Cross(x3_x1, x2_x1)) * 2;
            Vector3 grad1 = (grad2 * (-1.0f)) + (grad3 * (-1.0f));

            float gradientLengthSum = mySolver.invMass[x1[i]] * math.pow(grad1.magnitude, 2) + mySolver.invMass[x2[i]] * math.pow(grad2.magnitude, 2) + mySolver.invMass[x3[i]] * math.pow(grad3.magnitude, 2);

            //without damping
            float deltaLagrangeMultiplierNumerator1 = areaConstraint + (alpha * lambdas[i]);
            float deltaLagrangeMultiplierDenominator1 = gradientLengthSum + alpha;
            float deltaLambda = -deltaLagrangeMultiplierNumerator1 / deltaLagrangeMultiplierDenominator1;

            lambdas[i] += deltaLambda;

            //Debug.Log($"{grad1} / {grad2} / {grad3} * {mySolver.invMass[x1[i]]} * {deltaLambda}");
            mySolver.pos[x1[i]] += grad1 * mySolver.invMass[x1[i]] * deltaLambda / 20;
            mySolver.pos[x2[i]] += grad2 * mySolver.invMass[x2[i]] * deltaLambda / 20;
            mySolver.pos[x3[i]] += grad3 * mySolver.invMass[x3[i]] * deltaLambda / 20;
        }
    }
}



//public class CollisionConstraint : Constraint
//{
//    float restDistance = 0.06f;

//    float[] lambdas;
//    public CollisionConstraint(ClothXPBDSolver solver) : base(solver)
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
//        var solver = (ClothXPBDSolver)mySolver;
//        if (solver.collisions.Count <= 0)
//        {
//            return;
//        }
//        float alpha = stiff / (math.pow(dt, 2));


//        foreach (var pair in solver.collisions)
//        {
//            int index = pair.Key;
//            Vector3 distance = pair.Value;

//            float l = distance.magnitude;
//            float l_rest = restDistance;

//            float C = l - l_rest;
//            if (C >= 0.06)
//            {
//                continue;
//            }
//            //Debug.Log("wawa" +C.ToString());

//            //(xo-x1) * (1/|x0-x1|) = gradC
//            Vector3 gradC = distance.normalized;

  
//            float wTot = 1;

//            //lambda because |grad_Cn|^2 = 1 because if we move a particle 1 unit, the distance between the particles also grows with 1 unit, and w = w0 + w1
//            float deltalambda = -(C) / (wTot);
//            lambdas[index] += deltalambda;
//            //Move the vertices x = x + deltaX where deltaX = lambda * w * gradC
//            //Debug.Log($"{index}增加的距离为{deltalambda}");
//            mySolver.pos[index] += deltalambda * gradC;
//        }
//        //solver.collisions.Clear();
//    }
//}