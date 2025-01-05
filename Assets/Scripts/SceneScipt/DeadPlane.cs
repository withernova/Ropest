using System;
using UnityEngine;
using UnityEngine.Events;
using UnityEngine.Serialization;

public class DeadPlane : InteractiveBase
{
    [Serializable]
    public class PassVoid : UnityEvent { }

    [SerializeField] private PassVoid playerDead;


    public override void OnStart()
    {
        playerDead.Invoke();
    }

    public override void OnUpdate(float deltaTime)
    {
    }
}
