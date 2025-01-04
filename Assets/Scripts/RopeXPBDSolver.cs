using NUnit.Framework;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Xml.Serialization;
using Unity.Mathematics;
using UnityEngine;
using UnityEngine.UIElements;
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
        m_constraints.Add(new BendingAndTwistingConstraint(this));
        m_constraints.Add(new EdgeConstraint(this));
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

        float m = 1f;
        ghostInvMass = Enumerable.Repeat(m, pointPos.Count() - 1).ToArray();
        pointInvMass = Enumerable.Repeat(m, pointPos.Count()).ToArray();
        for (int i = 0; i < pointPos.Count(); i++)
        {
            pointInvMass[i] = i + 1;
            if(i != pointPos.Count() - 1)
                ghostInvMass[i] = i + 1;
        }
        for (int i = 0; i < pointPos.Count(); i++)
            new PointData(m);

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
            //if (i == grabPoint)
            //{
            //    continue;
            //}
            //v = (x - xPrev) / dt
            
            
            float f = -0.01f * pointInvMass[i];
            Vector3 deltaPos_G = 1f / 2 * dt * dt * gravity;
            Vector3 deltaPos_f = 1f / 2 * 0 * dt * dt * (pointPos[i] - prevPos[i]).normalized;
            pointPos[i] += deltaPos_f + deltaPos_G;
            if (i == ctrlIndex)
                pointPos[i] += move;
            
            // 用n帧的预测位置和n-1帧的实际位置得到n帧的预测速度
            vel[i] = (pointPos[i] - prevPos[i]) * oneOverdt;
            
            direcN[i] = Vector3.zero;
            foreach (var collider in Physics.OverlapCapsule(prevPos[i],pointPos[i], 0.03f))
            {
                // TODO: 添加排除自己的碰撞体
                if (collider.isTrigger) continue;
                Vector3 closestPoint = Physics.ClosestPoint(prevPos[i], collider,
                    collider.transform.position, collider.transform.rotation);
               
                //if(i == 10)
                //    DebugPoint(closestPoint, Color.white, 0.01f, dt);
                
                Physics.Raycast(prevPos[i], closestPoint - prevPos[i], out RaycastHit info, LayerMask.GetMask("Ignore Raycast"));
                Vector3 d0 = info.normal;

                vel[i] += Mathf.Clamp(-Vector3.Dot(vel[i], d0), 0, Mathf.Infinity) * d0;
                
                // var colliderTransform = collider.transform;
                // Vector3[] dVectors = new[]
                // {
                //     colliderTransform.up, -colliderTransform.up,
                //     colliderTransform.forward, -colliderTransform.forward,
                //     colliderTransform.right, -colliderTransform.right
                // };
                // float maxDotProduct = -1f; // 因为点乘值可能是负数，初始化为-1
                // Vector3 closestVector = Vector3.zero;
                // foreach (Vector3 di in dVectors)
                // {
                //     float dotProduct = Vector3.Dot(d0.normalized, di.normalized); // 计算点乘（归一化以比较方向）
                //
                //     // 找到点乘值最接近1的向量
                //     if (dotProduct > maxDotProduct)
                //     {
                //         maxDotProduct = dotProduct;
                //         closestVector = di;
                //     }
                // }
                // d0 = closestVector;
                
                // if(i == 10)
                //     Debug.Log("处理后法向量:"+ d0);
                
                //direcN[i] += d0;
            }
            pointPos[i] = prevPos[i] + vel[i] * dt;

            // 往前对于pointPos的更新是改变n帧的预测位置
            
            if (i < pointPos.Count() - 1)
                ghostVels[i] = (ghostPos[i] - ghostPrev[i]) * oneOverdt;

            Enforce();

        }
    }



     protected override void PreSolve(float dt, Vector3 g)
    {
        forces.ForEach((f) => { vel[f.Key] += dt * f.Value; });
        forces.Clear();
        for (int i = 0; i < pointPos.Count(); i++)
        {
            if (pointInvMass[i] == 0f)
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

            
            // n+1帧：往前的pointpos是第n帧的实际位置
            
            prevPos[i] = pointPos[i];
            pointPos[i] += vel[i] * dt;
            
            //pointPos[i] = Vector3.Lerp(pointPos[i], pointPos[i] + vel[i] * dt, 0.8f);

            if (i < pointPos.Count() - 1)
            {

                ghostPrev[i] = ghostPos[i];
                ghostPos[i] += ghostVels[i] * dt;
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

    
    private Vector3 enforceMove = new Vector3();
   
    public void Enforce()
    {
        if (enforceMove.magnitude > 1e-6f)
        {
            pointPos[ctrlIndex] += enforceMove;
            enforceMove = new Vector3();
        }
        for(int i=0;i<PointData.datas.Count;i++)
        {
            var data = PointData.datas[i]; 
            if (data.interactiveItem is InteractiveGrab)
                //((InteractiveGrab)data.interactiveItem).SetV(vel[i]) ;
                data.interactiveItem.transform.position = pointPos[ctrlIndex];

        }
    }
    
    public void Swing(InteractiveSwing target)
    {
        Debug.Log("swing");
        enforceMove = target.GetTargetPos(pointPos[ctrlIndex]);
        PointData.datas[ctrlIndex].SetActive(target);
        pointInvMass[ctrlIndex] = 0;
    }

    public void Grab(InteractiveGrab target)
    {
        Debug.Log("grab");
        PointData.datas[ctrlIndex].SetActive(target);

        //先设置更新的位移
        //enforceMove = target.transform.position - pointPos[ctrlIndex];
        target.transform.position = pointPos[ctrlIndex];

        //然后设置质量
        pointInvMass[ctrlIndex] += target.GetMass(); 
    }

    public void LoseControl()
    {
        pointInvMass[ctrlIndex] = PointData.datas[ctrlIndex].Reset();
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

        for (int e = 0; e < ghostPos.Count() - 1; e++)
        {
            Vector3 v1 = pointPos[e];
            Vector3 v2 = pointPos[e + 1];
            Vector3 v3 = pointPos[e + 2];

            Vector3 vm = 0.5f * (v1 + v2);
            Vector3 vml = 0.5f * (v2 + v3);

            Vector3 d3f = (v2 - v1).normalized;
            Vector3 d2f = Vector3.Cross(d3f,(ghostPos[e] - v1).normalized).normalized;
            Vector3 d1f = Vector3.Cross(d2f, d3f).normalized;

            Vector3 d3l = (v3 - v2).normalized;
            Vector3 d2l = Vector3.Cross(d3l,(ghostPos[e + 1] - v2).normalized).normalized;
            Vector3 d1l = Vector3.Cross(d2l, d3l).normalized;   

            Matrix4x4 De1 = Matrix4x4.identity;
            Matrix4x4 De2 = Matrix4x4.identity;
            De1.SetColumn(0, d1f);
            De1.SetColumn(1, d2f);
            De1.SetColumn(2, d3f);
            De2.SetColumn(0, d1l);
            De2.SetColumn(1, d2l);
            De2.SetColumn(2, d3l);

            Matrix4x4 rotation = De2*Matrix4x4.Transpose(De1) ;
            float trace = rotation.m00 + rotation.m11 + rotation.m22;
            float theta = Mathf.Acos((trace - 1) / 2f);
            Vector3 n = (theta > 1e-6f) ? new Vector3(
                (rotation.m21 - rotation.m12) / (2f * Mathf.Sin(theta)),
                (rotation.m02 - rotation.m20) / (2f * Mathf.Sin(theta)),
                (rotation.m10 - rotation.m01) / (2f * Mathf.Sin(theta))
            ) : Vector3.zero;

            float le = (vml - vm).magnitude;
            int interplor = 10;
            List<Vector3> sPosition = new List<Vector3>();
            List<Matrix4x4> interpolatedFrames = new List<Matrix4x4>();
            float segmentLength = le / (interplor-1);
            Vector3 currentPosition = vm;
            sPosition.Add(currentPosition);


            for (int i = 0; i < interplor; i++)
            {
                float r = ((i * segmentLength)) / le;
                Vector3 scaledTheta = theta * n.normalized * r;
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

                Matrix4x4 interpolatedFrame = rotationInterpolated * De1;
                interpolatedFrames.Add(interpolatedFrame);

                if (i > 0)
                {
                    Vector3 d3 = (Vector3)interpolatedFrame.GetColumn(2).normalized;
                    Vector3 d2 = (Vector3)interpolatedFrame.GetColumn(1).normalized;
                    Vector3 d1 = (Vector3)interpolatedFrame.GetColumn(0).normalized;

                    //Debug.Log($"Segment {e}, Interpolated {i}: d1={d1}, d2={d2}, d3={d3}");
                    currentPosition = currentPosition + d3 * segmentLength;
                    sPosition.Add(currentPosition);
                }
            }
                //sPosition.Add(vml);
            // -----------------------------------

            //string log = "";
            //Debug.Log("_________________________________");
            //Debug.Log($"非插值前d3{d3f} ");
            //interpolatedFrames.ForEach(frame => log += $"插值d3:{(Vector3)frame.GetColumn(2)} ");
            //Debug.Log(log);
            //Debug.Log($"非插值后d3{d3l}");

            //Debug.Log("_________________________________");

            //for (int j = 0; j < subdivision; ++j)
            //{
            //    float angle = 2 * math.PI * j / subdivision;
            //    Vector3 vertex = vm + radius * (math.cos(angle) * (Vector3)d1f.normalized + math.sin(angle) * (Vector3)d2f.normalized);
            //    newpos.Add(vertex);
            //}

            // 将截面点存储到 allSections 中
            List<Vector3> sectionPoints = new List<Vector3>();
            for (int i = 0; i < interplor-1; i++)
            {
                for (int j = 0; j < subdivision; ++j)
                {

                    float angle = 2 * math.PI * j / subdivision;
                    Vector3 vertex = sPosition[i] + radius * (math.cos(angle) * (Vector3)interpolatedFrames[i].GetColumn(0).normalized + math.sin(angle) * (Vector3)interpolatedFrames[i].GetColumn(1).normalized);
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