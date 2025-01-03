using System.Collections.Generic;
using System.Linq;
using UnityEngine;

public class SolverInitData
{

}

public abstract class XPBDSolver
{
    protected List<Constraint> m_constraints = new List<Constraint>();
    protected int numSubSteps = 7;
    public Vector3 gravity;

    public Vector3[] pos;
    public Vector3[] prevPos, vel;
    public Vector3[] direcN;

    public float[] invMass;
    public Triangle[] triangles;


    public int numParticles;

    public List<KeyValuePair<int, Vector3>> forces = new List<KeyValuePair<int, Vector3>>();
    public int grabPoint = -1;

    protected float grabInvMass;

    public void InitByMesh(Mesh mesh, Vector3 g, SolverInitData data)
    {
        gravity = g;
        pos = mesh.vertices;
        prevPos = mesh.vertices;
        numParticles = pos.Length;
        vel = new Vector3[numParticles];
        direcN = Enumerable.Repeat(Vector3.zero, pos.Length).ToArray();
        

        List<Triangle> triangleList = new List<Triangle>();
        //float totalMass = 1;
        //float totalS = 0;
        invMass = Enumerable.Repeat(0f, pos.Length).ToArray();
        for (int i = 0; i < mesh.triangles.Length;)
        {
            int i0 = i;
            Vector3[] points = new Vector3[3];
            int[] pointsId = new int[3];
            for (; i < i0 + 3; i++)
            {
                points[i - i0] = mesh.vertices[mesh.triangles[i]];
                pointsId[i - i0] = mesh.triangles[i];
            }

            float a = (points[0] - points[1]).magnitude;
            float b = (points[2] - points[1]).magnitude;
            float c = (points[0] - points[2]).magnitude;
            float p = (a + b + c) / 2;
            float s = Mathf.Sqrt(p * (p - a) * (p - b) * (p - c));

            triangleList.Add(new Triangle(pointsId));
            //totalS += s;
            invMass[pointsId[0]] += s / 3;
            invMass[pointsId[1]] += s / 3;
            invMass[pointsId[2]] += s / 3;
        }

        triangles = triangleList.ToArray();
        for (int i = 0; i < invMass.Length; i++)
        {
            invMass[i] = 1f / invMass[i];
        }

        InitSolver(mesh, data);

        AddConstraints();
    }

    public void AddForce(int pointIndex, Vector3 force)
    {
        forces.Add(new KeyValuePair<int, Vector3>(pointIndex, force));
    }

    public void Simulate(float dt)
    {
        float sdt = dt;

        m_constraints.ForEach((c) => c.ResetLambda());

        PreSolve(dt, gravity);
        for (int step = 0; step < numSubSteps; step++)
        {

            SolveConstraints(dt);

        }
        PostSolve(dt);
    }

    // 处理重力约束
    protected virtual void PreSolve(float dt, Vector3 g)
    {
        g = Vector3.zero;
        forces.ForEach((f) => { vel[f.Key] += dt * f.Value; Debug.Log(f); });
        forces.Clear();
        for (int i = 0; i < numParticles; i++)
        {
            if (invMass[i] == 0f)
            {
                continue;
            }

            prevPos[i] = pos[i];
            //pos[i] = Vector3.Lerp(pos[i], pos[i] + vel[i] * dt, 1f);
            pos[i] += vel[i] * dt;
        }
    }

    // 处理其他约束
    private void SolveConstraints(float dt)
    {
        foreach (Constraint constraint in m_constraints)
        {
            constraint.SolveConstraint(dt);
        }
    }

    // 更新速度
    protected virtual void PostSolve(float dt)
    {
        float oneOverdt = 1f / dt;

        //For each particle
        for (int i = 0; i < numParticles; i++)
        {
            if (invMass[i] == 0f)
            {
                continue;
            }
            //v = (x - xPrev) / dt
            float f = -0.01f * invMass[i];
            float g = -9.8f;
            Vector3 deltaPos_G = 1f / 2 * g * dt * dt * new Vector3(0, 1, 0);
            Vector3 deltaPos_f = 1f / 2 * f * dt * dt * (pos[i] - prevPos[i]).normalized;
            pos[i] += deltaPos_f + deltaPos_G;
            // if (i > numParticles / 2)
            // {
            //     pos[i] += 1f / 2 * a * dt * dt * new Vector3(0, 1, 0);
            // }
            vel[i] = (pos[i] - prevPos[i]) * oneOverdt;
        }
    }

    public virtual void StartGrab(Vector3 grabPos, Transform trans)
    {
        grabPoint = FindClosestPoint(grabPos, trans);

        grabInvMass = invMass[grabPoint];
        invMass[grabPoint] = 0;
        vel[grabPoint] = Vector3.zero;

        Debug.Log(grabPoint);
    }

    public virtual void OnGrabbing(Vector3 grabPos, Transform trans)
    {
        //Debug.DrawLine(grabPos, trans.TransformPoint(pos[grabPoint]));
        pos[grabPoint] = trans.InverseTransformPoint(grabPos);
    }

    public virtual void EndGrab(Vector3 grabPos, Transform trans)
    {
        invMass[grabPoint] = grabInvMass;
        grabPoint = -1;
    }

    public virtual int FindClosestPoint(Vector3 toFindPos, Transform trans)
    {
        float min = pos.ToList().Min(p => (toFindPos - trans.TransformPoint(p)).magnitude);
        return pos.ToList().FindIndex(p => (toFindPos - trans.TransformPoint(p)).magnitude == min);
    }
    protected abstract void AddConstraints();
    protected abstract void InitSolver(Mesh mesh, SolverInitData data);
}