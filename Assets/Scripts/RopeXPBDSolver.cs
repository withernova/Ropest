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

    public Quaternion[] pointQ;

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
        //m_constraints.Add(new Volume2Constraint(this));
        m_constraints.Add(new BendingAndTwistingConstraint(this));
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

        pointQ = new Quaternion[pointPos.Count() - 1];
        Vector3 from = new Vector3(0, 0, 1);
        for (int i = 0; i < pointPos.Count() - 1; i++)
        {
            Vector3 to = (points[i + 1] - points[i]).normalized;
            Quaternion dq = Quaternion.FromToRotation(from, to);
            if (i == 0) pointQ[i] = dq;
            else pointQ[i] = dq * pointQ[i - 1];
            from = to;
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

            //if (i == 0 || i == pointPos.Count() - 1)
            //{
            //    pointPos[i] = prevPos[i];
            //    continue;
            //}
            //vel[i] += g * dt;

            //prevPos[i] = pointPos[i];
            //pointPos[i] = Vector3.Lerp(pointPos[i], pointPos[i] + vel[i] * dt, 0.8f);

            //if (i < pointPos.Count() - 1)
            //{

            //    ghostPrev[i] = ghostPos[i];
            //    ghostPos[i] = Vector3.Lerp(ghostPos[i], ghostPos[i] + ghostVels[i] * dt, 0.8f);
            //}

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

        for (int e = 0; e < ghostPos.Count(); e++)
        {
            Vector3 v1 = pointPos[e];
            Vector3 v2 = pointPos[e + 1];
            //Quaternion q = pointQ[e];

            Vector3 vm = 0.5f * (v1 + v2);

            //Vector3 d1 = q * new Vector3(1, 0, 0); //* scale;
            //Vector3 d2 = q * new Vector3(0, 1, 0);//* scale;
            //Vector3 d3 = q * new Vector3(0, 0, 1); //* scale;
            Vector3 d2 = (ghostPos[e] - vm).normalized;
            Vector3 d1 = Vector3.Cross(vm.normalized, d2).normalized;

            //Debug.Log($"{d1} {d2} {d3}");
            for (int j = 0; j < subdivision; ++j)
            {
                float angle = 2 * math.PI * j / subdivision;
                Vector3 vertex = v1 + radius * (math.cos(angle) * d1 + math.sin(angle) * d2);
                newpos.Add(vertex);
            }

            for (int i = 0; i < 10; i++)
            {
                if (e != pointPos.Count() - 2)
                {
                    Vector3 nd2 = (ghostPos[e + 1] - 0.5f * (pointPos[e + 1] + pointPos[e + 2])).normalized;
                    Vector3 nd1 = Vector3.Cross(0.5f * (pointPos[e + 1] + pointPos[e + 2]), d2).normalized;

                    for (int j = 0; j < subdivision; ++j)
                    {
                        float angle = 2 * math.PI * j / subdivision;
                        Vector3 vertex = Vector3.Lerp(v1, v2, i * 1f/10f) + radius * (math.cos(angle) * Vector3.Lerp(d1, nd1, i * 1f/10f) + math.sin(angle) * Vector3.Lerp(d2, nd2, i * 1f / 10f));
                        newpos.Add(vertex);
                    }
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

            if (e == pointPos.Count() - 2)
            {
                for (int j = 0; j < subdivision; ++j)
                {
                    float angle = 2 * math.PI * j / subdivision;
                    Vector3 vertex = pointPos.Last() + radius * (math.cos(angle) * d1 + math.sin(angle) * d2);
                    newpos.Add(vertex);
                }
            }
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