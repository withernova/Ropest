using System;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;
public class ClothSolverInitData : SolverInitData
{
    public int subdivision;

    public ClothSolverInitData(int subdivision)
    {
        this.subdivision = subdivision;
    }
}

public class ClothXPBDSolver : XPBDSolver
{
    public int[][] sections;
    public int subdivision;

    public float sumVolume = 0f;

    protected override void AddConstraints()
    {
        m_constraints.Add(new DistanceConstraint(this));
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
            pos[i] += vel[i] * dt;
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
            float f = -0.01f * invMass[i];
            Vector3 deltaPos_G = 1f / 2 * dt * dt * gravity;
            Vector3 deltaPos_f = 1f / 2 * f * dt * dt * (pos[i] - prevPos[i]).normalized;
            pos[i] += deltaPos_f + deltaPos_G;
            vel[i] = (pos[i] - prevPos[i]) * oneOverdt;
            
            direcN[i] = direcN[i].normalized;
            
            //Debug.Log(direcN);
            vel[i] += Mathf.Clamp(-Vector3.Dot(vel[i], direcN[i]), 0, Mathf.Infinity)  * direcN[i];
            pos[i] = prevPos[i] + vel[i] * dt;
        }
    }
}
