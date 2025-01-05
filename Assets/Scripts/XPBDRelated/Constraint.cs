using System.Collections.Generic;
using System.Linq;
using TMPro;
using Unity.Mathematics;
using Unity.VisualScripting;
using UnityEngine;
using static UnityEditor.PlayerSettings;
using static UnityEditor.ShaderGraph.Internal.KeywordDependentCollection;
using static ClothXPBDSolver;


public abstract class Constraint
{
    public XPBDSolver mySolver;

    public float stiff = 0.6f;

    public Constraint(XPBDSolver solver)
    {
        mySolver = solver;
    }

    public abstract void SolveConstraint(float dt);

    public abstract void ResetLambda();
}