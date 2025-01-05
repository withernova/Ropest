using System;
using UnityEngine;
using UnityEngine.Events;

public class TipsArea : InteractiveBase
{
    [Serializable]
    public class PassVoid : UnityEvent { }

    [SerializeField] private DeadPlane.PassVoid playerEnter;
    
    public override void OnStart()
    {
        
    }

    public override void OnUpdate(float deltaTime)
    {
        
    }
}
