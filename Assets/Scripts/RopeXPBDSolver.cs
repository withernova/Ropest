using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using Unity.Mathematics;
using UnityEngine;
using static UnityEditor.Searcher.SearcherWindow.Alignment;

public class RopeSolverInitData : SolverInitData
{
    public int subdivision;
    public int segments;
    public RopeSolverInitData(int subdivision, int segments)
    {
        this.subdivision = subdivision;
        this.segments = segments;
    }
}
public class RopeXPBDSolver : XPBDSolver, IControllable
{
    public Vector3[] pointPos;
    public Vector3[] ghostPos;
    public Vector3[] ghostPrev;
    public Vector3[] ghostVels;


    public float[] ghostInvMass;
    public float[] pointInvMass;
    public int[][] sections;
    public int segments;
    public float radius = 0.01f;
    public int subdivision;
    public float[] length;

    public float ghostDistance = 1f;

    public int ctrlIndex = 0;

    private Vector3 move;
    
    public void SetMove(Vector3 move)
    {
        this.move = move;
    }

    protected override void AddConstraints()
    {
        m_constraints.Add(new EdgeConstraint(this));
        //m_constraints.Add(new DoubleDistanceConstraint(this));

        m_constraints.Add(new BendingAndTwistingConstraint(this));
        //m_constraints.Add(new TargetConstraint(this));
        //m_constraints.Add(new Volume2Constraint(this));
        //m_constraints.Add(new BendTwistConstraint(this));
    }

    protected override void InitSolver(Mesh mesh, SolverInitData data)
    {
        numSubSteps = 10;
        RopeSolverInitData ropeData = data as RopeSolverInitData;
        subdivision = ropeData.subdivision;

        var sectionList = new List<int[]>();
        var points = new List<Vector3>();

        int numSections = (mesh.vertexCount - 2) / subdivision; // 计算截面数量
        for (int i = 0; i < numSections; i++)
        {
            List<int> sectionPoints = new List<int>();
            Vector3 center = new Vector3();
            for (int j = 0; j < subdivision; j++)
            {
                sectionPoints.Add(i * subdivision + 1 + j);
                center += pos[i * subdivision + 1 + j];
            }
            sectionList.Add(sectionPoints.ToArray());
            center /= subdivision;
            points.Add(center);
        }
        sections = sectionList.ToArray();
        pointPos = points.ToArray();
        prevPos = points.ToArray();
        ghostPos = new Vector3[pointPos.Count() - 1];

        ghostVels = new Vector3[pointPos.Count() - 1];

        for (int i = 0; i < numSections - 1; i++)
        {
            ghostPos[i] = (pointPos[i] + pointPos[i + 1]) / 2;
            ghostPos[i].y += ghostDistance;
        }
        ghostPrev = ghostPos;

        ghostInvMass = Enumerable.Repeat(100f, pointPos.Count() - 1).ToArray();
        pointInvMass = Enumerable.Repeat(100f, pointPos.Count()).ToArray();

        length = new float[ghostPos.Count()];
        for (int i = 0; i < ghostPos.Count(); i++)
        {
            length[i] = (pointPos[i] - pointPos[i + 1]).magnitude;
        }


    }

    protected override void PostSolve(float dt)
    {
        float oneOverdt = 1f / dt;

        //For each particle
        for (int i = 0; i < pointPos.Count(); i++)
        {
            if (i == grabPoint)
            {
                continue;
            }
            //v = (x - xPrev) / dt
            float f = -0.0001f * 1/invMass[i];
            float g = -9.8f;
            Vector3 deltaPos_G = 1f / 2 * g * dt * dt * Vector3.up;
            Vector3 deltaPos_f = 1f / 2 * f * dt * dt * (pointPos[i] - prevPos[i]).normalized;
            pointPos[i] += deltaPos_f +  deltaPos_G;
            vel[i] = (pointPos[i] - prevPos[i]) * oneOverdt;
            vel[i] += Mathf.Clamp(-Vector3.Dot(vel[i], direcN[i]), 0, Mathf.Infinity) * direcN[i];
            pointPos[i] = prevPos[i] + vel[i] * dt;
            if (i < pointPos.Count() - 1)
                ghostVels[i] = (ghostPos[i] - ghostPrev[i]) * oneOverdt;
        }
    }

    protected override void PreSolve(float dt, Vector3 g)
    {
        forces.ForEach((f) => { vel[f.Key] += dt * f.Value; });
        forces.Clear();
        for (int i = 0; i < pointPos.Count(); i++)
        {
            if (invMass[i] == 0f)
            {
                prevPos[i] = pointPos[i];
                continue;
            }

            if (i == 0 || i == pointPos.Count() - 1)
            {
                pointPos[i] = prevPos[i];
                continue;
            }
            //vel[i] += g * dt;

            prevPos[i] = pointPos[i];
            pointPos[i] = Vector3.Lerp(pointPos[i], pointPos[i] + vel[i] * dt, 0.8f);

            if (i < pointPos.Count() - 1)
            {

                ghostPrev[i] = ghostPos[i];
                ghostPos[i] = Vector3.Lerp(ghostPos[i], ghostPos[i] + ghostVels[i] * dt, 0.8f);
            }

        }
    }

    public override int FindClosestPoint(Vector3 toFindPos, Transform trans)
    {
        float min = pointPos.ToList().Min(p => (toFindPos - trans.TransformPoint(p)).magnitude);
        return pointPos.ToList().FindIndex(p => (toFindPos - trans.TransformPoint(p)).magnitude == min);
    }

    public override void OnGrabbing(Vector3 grabPos, Transform trans)
    {
        //Debug.DrawLine(grabPos, trans.TransformPoint(pos[grabPoint]));
        for (int i = 0; i < ghostPos.Count(); i++)
        {
            Debug.DrawLine(trans.TransformPoint(ghostPos[i]), trans.TransformPoint(ghostPos[i]) + Vector3.one * 0.1f);
        }
        pointPos[grabPoint] = trans.InverseTransformPoint(grabPos);
    }

    //public override void StartGrab(Vector3 grabPos, Transform trans)
    //{
    //    grabPoint = FindClosestPoint(grabPos, trans);

    //    grabInvMass = invMass[grabPoint];
    //    vel[grabPoint] = Vector3.zero;

    //    Debug.Log(grabPoint);
    //}

    //public override void EndGrab(Vector3 grabPos, Transform trans)
    //{
    //    invMass[grabPoint] = grabInvMass;
    //    grabPoint = -1;
    //}

    public void TryRendering(Mesh mesh)
    {
        //更新pos
        List<Vector3> newpos = new List<Vector3>();
        newpos.Add(pointPos[0]);

        for (int e = 0; e < ghostPos.Count()-1; e++)
        {
            Vector3 v1 = pointPos[e];
            Vector3 v2 = pointPos[e + 1];
            Vector3 v3 = pointPos[e + 2];

            //Quaternion q = pointQ[e];

            Vector3 vm = 0.5f * (v1 + v2);
            Vector3 vml = 0.5f * (v2 + v3);

            //Vector3 d1 = q * new Vector3(1, 0, 0); //* scale;
            //Vector3 d2 = q * new Vector3(0, 1, 0);//* scale;
            //Vector3 d3 = q * new Vector3(0, 0, 1); //* scale;
            Vector3 d3f = (v2 - v1).normalized;
            Vector3 d2f = Vector3.Cross((ghostPos[e] - v1),d3f).normalized;
            Vector3 d1f = Vector3.Cross(d3f, d2f).normalized;
            Matrix4x4 De1 = new Matrix4x4();
            Vector3 d3l = (v3 - v2).normalized;
            Vector3 d2l = Vector3.Cross((ghostPos[e+1] - v2),d3l).normalized;
            Vector3 d1l = Vector3.Cross(d3l, d2l).normalized;
            Matrix4x4 De2 = new Matrix4x4();

            De1.SetColumn(0,d1f);
            De1.SetColumn(1,d2f);
            De1.SetColumn(2,d3f);
            De2.SetColumn(0, d1l);
            De2.SetColumn(1, d2l);
            De2.SetColumn(2, d3l);

            Matrix4x4 rotation = Matrix4x4.Transpose(De1) * De2;

            float trace = rotation.m00 + rotation.m11 + rotation.m22;
            float theta = Mathf.Acos((trace - 1) / 2f);
            Vector3 n;

            if (theta > 1e-6f)  // 如果 θ ≠ 0
            {
                float factor = 1f / (2f * Mathf.Sin(theta));
                n = new Vector3(
                    (rotation.m21 - rotation.m12) * factor,
                    (rotation.m02 - rotation.m20) * factor,
                    (rotation.m10 - rotation.m01) * factor
                );
            }
            else
            {
                n = Vector3.zero;  // 如果 θ = 0，旋转轴无意义
            }

            float le = (vml - vm).magnitude;
            int interplor = 2;

            List<Vector3> sPostion = new List<Vector3>();
            List<float> slist= new List<float>();
            for(int i = 0; i < interplor; i++)
            {
                float t = i / (float)(interplor - 1);

                // 线性插值计算位置
                sPostion.Add(Vector3.Lerp(vm, vml, t));
                slist.Add((i+1)* le / interplor);

            }

            //下面的罗德里格斯公式可能有问题
            List<Matrix4x4> interpolatedFrames = new List<Matrix4x4>();
            int enumNum = interplor;
            float segmentLength = le / (interplor - 1);
            for (int i = 0; i < enumNum; i++)
            {
                // 限制弧长范围在 [0, le]
                float s = Mathf.Clamp(i * segmentLength, 0, le);
                float r = (s - (le / 2f)) / le;

                // 确保比例 r 在 [0, 1] 范围内
                //if (r < 0 || r > 1)
                //{
                //    interplor--;
                //}

                Vector3 scaledTheta = theta * n * r;

                // 使用 Rodrigues 公式计算旋转矩阵
                float cosR = Mathf.Cos(scaledTheta.magnitude);
                float sinR = Mathf.Sin(scaledTheta.magnitude);
                Vector3 axis = scaledTheta.normalized;

                Matrix4x4 rotationInterpolated = Matrix4x4.identity;
                rotationInterpolated.m00 = cosR + axis.x * axis.x * (1 - cosR);
                rotationInterpolated.m01 = axis.x * axis.y * (1 - cosR) - axis.z * sinR;
                rotationInterpolated.m02 = axis.x * axis.z * (1 - cosR) + axis.y * sinR;

                rotationInterpolated.m10 = axis.y * axis.x * (1 - cosR) + axis.z * sinR;
                rotationInterpolated.m11 = cosR + axis.y * axis.y * (1 - cosR);
                rotationInterpolated.m12 = axis.y * axis.z * (1 - cosR) - axis.x * sinR;

                rotationInterpolated.m20 = axis.z * axis.x * (1 - cosR) - axis.y * sinR;
                rotationInterpolated.m21 = axis.z * axis.y * (1 - cosR) + axis.x * sinR;
                rotationInterpolated.m22 = cosR + axis.z * axis.z * (1 - cosR);

                // 插值材料帧
                Matrix4x4 interpolatedFrame = rotationInterpolated * De1;
                interpolatedFrames.Add(interpolatedFrame);
            }
            // -----------------------------------

            //string log = "";
            //Debug.Log("_________________________________");
            //Debug.Log($"非插值前d1{d1f} 非插值前d2{d2f}");
            //interpolatedFrames.ForEach(frame => log += $"插值d1:{(Vector3)frame.GetColumn(0)}插值d2:{(Vector3)frame.GetColumn(1)} ");
            //Debug.Log(log);
            //Debug.Log($"非插值后d1{d1l} 非插值后d2{d2l}");

            //Debug.Log("_________________________________");

            //for (int j = 0; j < subdivision; ++j)
            //{
            //    float angle = 2 * math.PI * j / subdivision;
            //    Vector3 vertex = vm + radius * (math.cos(angle) * (Vector3)d1f.normalized + math.sin(angle) * (Vector3)d2f.normalized);
            //    newpos.Add(vertex);
            //}

            // 将截面点存储到 allSections 中
            List<Vector3> sectionPoints = new List<Vector3>();
            for (int i = 1; i < interplor; i++)
            {
                for (int j = 0; j < subdivision; ++j)
                {
                    float angle = 2 * math.PI * j / subdivision;
                    Vector3 vertex = sPostion[i] + radius * (math.cos(angle) * (Vector3)interpolatedFrames[i].GetColumn(0).normalized + math.sin(angle) * (Vector3)interpolatedFrames[i].GetColumn(1).normalized);
                    sectionPoints.Add(vertex);
                    newpos.Add(vertex);

                }
            }

            //for (int j = 0; j < subdivision; ++j)
            //{
            //    float angle = 2 * math.PI * j / subdivision;
            //    Vector3 vertex = v2 + radius * (math.cos(angle) * d1 + math.sin(angle) * d2);
            //    newpos.Add(vertex);
            //}

            //for (int j = 0; j < subdivision; ++j)
            //{
            //    float angle = 2 * math.PI * j / subdivision;
            //    Vector3 vertex = v2 + radius * (math.cos(angle) * d1 + math.sin(angle) * d2);
            //    newpos.Add(vertex);
            //}

            //if (e == pointPos.Count() - 2)
            //{
            //    for (int j = 0; j < subdivision; ++j)
            //    {
            //        float angle = 2 * math.PI * j / subdivision;
            //        Vector3 vertex = pointPos.Last() + radius * (math.cos(angle) * d1 + math.sin(angle) * d2);
            //        newpos.Add(vertex);
            //    }
            //}
        }

        newpos.Add(pointPos.Last());

        //pos = newpos.ToArray();
        // 生成网格
        GenerateMeshCosserat(newpos, mesh, subdivision);
    }

    void GenerateMeshCosserat(List<Vector3> vers, Mesh mesh, int subdivisions)
    {
        int numPoints = (vers.Count() - 2) / subdivisions;
        int sub = subdivisions;

        List<Vector3> vertices = vers;
        List<int> triangles = new List<int>();
        List<Vector2> uvs = new List<Vector2>();

        int seg = numPoints - 1;

        // 索引
        for (int i = 0; i < seg; i++)
        {
            if (i == 0)
            {
                for (int j = 0; j < sub; j++)
                {
                    int nextJ = (j + 1) % sub;
                    triangles.Add(0);
                    triangles.Add(1 + nextJ);
                    triangles.Add(1 + j);
                }
            }
            for (int j = 0; j < sub; j++)
            {
                int nextJ = (j + 1) % sub;
                triangles.Add(i * sub + 1 + j);
                triangles.Add(i * sub + 1 + nextJ);
                triangles.Add((i + 1) * sub + 1 + j);

                triangles.Add((i + 1) * sub + 1 + j);
                triangles.Add(i * sub + 1 + nextJ);
                triangles.Add((i + 1) * sub + 1 + nextJ);
            }
        }
        for (int j = 0; j < sub; j++)
        {
            int nextJ = (j + 1) % sub;
            triangles.Add(seg * sub + 1 + j);
            triangles.Add(seg * sub + 1 + nextJ);
            triangles.Add((seg + 1) * sub + 1);
        }


        // uv
        uvs.Add(new Vector2(0, 0));
        for (int i = 0; i <= seg; i++)
        {
            for (int j = 0; j < sub; j++)
            {
                float u, v;
                if (j <= sub / 2)
                {
                    u = j / sub * 2;
                }
                else
                {
                    u = (sub - j) / sub * 2;
                }
                v = i;
                uvs.Add(new Vector2(u, v));
            }
        }
        uvs.Add(new Vector2(0, 0));

        mesh.Clear();
        mesh.SetVertices(vertices);
        mesh.SetTriangles(triangles, 0);
        mesh.uv = uvs.ToArray();
        mesh.RecalculateNormals();
        mesh.RecalculateBounds();
    }
}