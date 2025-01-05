using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Cryptography;
using UnityEngine;
public class ClothSolverInitData : SolverInitData
{
    public int subdivision;
    public SphereCollider cc;
    public ClothSolverInitData(int subdivision, SphereCollider cc)
    {
        this.subdivision = subdivision;
        this.cc = cc;
    }
}

public class ClothXPBDSolver : XPBDSolver
{
    public int[][] sections;
    public int subdivision;

    public float sumVolume = 0f;
    private SphereCollider cc;

    protected override void AddConstraints()
    {
        m_constraints.Add(new DistanceConstraint(this));
        //m_constraints.Add(new CollisionConstraint(this));
        //m_constraints.Add(new VolumeConstraint(this));
        //m_constraints.Add(new SectionVolume2Constraint(this));
        //m_constraints.Add(new SectionAreaConstriant(this));
        //m_constraints.Add(new SectionDistanceConstraint(this));
        //m_constraints.Add(new RoundConstraint(this));
        //m_constraints.Add(new ElasticConstraint(this));
        //m_constraints.Add(new BendingConstraint(this));
    }

    protected override void InitSolver(Mesh mesh, SolverInitData data)
    {
        ClothSolverInitData clothData = data as ClothSolverInitData;
        subdivision = clothData.subdivision;
        cc = clothData.cc;
        cc.radius = 0.02f;

        var sectionList = new List<int[]>();

        int numSections = (mesh.vertexCount - 2) / subdivision; // 计算截面数量
        for (int i = 0; i < numSections; i++)
        {
            List<int> sectionPoints = new List<int>();
            for (int j = 0; j < subdivision; j++)
            {
                sectionPoints.Add(i * subdivision + 1 + j);
            }
            sectionList.Add(sectionPoints.ToArray());
        }
        sections = sectionList.ToArray();
    }

    protected override void PreSolve(float dt, Vector3 g)
    {
        forces.ForEach((f) => { vel[f.Key] += dt * f.Value; Debug.Log(f); });
        forces.Clear();
        for (int i = 0; i < numParticles; i++)
        {
            if (invMass[i] == 0f)
            {
                continue;
            }

            vel[i] += g * dt;

            prevPos[i] = pos[i];
            pos[i] += vel[i] * dt;
            //Debug.Log(vel[i]);
        }
    }

    protected override void PostSolve(float dt)
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

            foreach (var collider in Physics.OverlapSphere(pos[i], 0.03f, ~LayerMask.GetMask("Ignore Raycast")))
            {
                if (collider.isTrigger) continue;

                SetSphereTo(pos[i]);

                Vector3 otherPosition = collider.gameObject.transform.position;
                Quaternion otherRotation = collider.gameObject.transform.rotation;
                Vector3 direction;
                float distance;

                bool overlapped = Physics.ComputePenetration(
                    cc, cc.transform.position, cc.transform.rotation,
                    collider, otherPosition, otherRotation,
                    out direction, out distance
                );
                Vector3 correctionVector = distance * direction;
                //the most inelegant way, but time is short
                if (overlapped)
                {
                    pos[i] += correctionVector;
                    if (collider.name == "collider") Debug.Log(correctionVector);
                    FrictionForNative(i, distance, direction);
                }
            }
            vel[i] = (pos[i] - prevPos[i]) * oneOverdt;
        }
    }

    void SetSphereTo(Vector3 start)
    {
        var capParent = cc.transform;

        //Set mid of capsule to mid of in-between vector
        capParent.position = start;
    }


    void FrictionForNative(int i, float distance, Vector3 direction)
    {
        ref var p1 = ref pos[i];
        ref var p1_ = ref prevPos[i];

        var p1_p1 = p1 - p1_;

        var t1 = Vector3.ProjectOnPlane(p1_p1, direction);

        if (t1.magnitude < 0.01f)
        {
            p1_ = p1;
            return;
        }

        p1_ += p1_p1.normalized * 0.25f * distance;
    }

}
