using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Events;

public abstract class InteractiveBase : MonoBehaviour
{
    public float lastingTime;
    public bool activating;
    public float radius;
    public UnityAction OnInteractive;
    public abstract void OnStart();
    public abstract void OnUpdate(float deltaTime);

    public virtual void EndInteractive()
    {
       
    }

    private void OnDrawGizmos()
    {
        Gizmos.color = new Color(0, 1, 0, 0.2f);
        Gizmos.DrawSphere(transform.position, radius);
    }

}
