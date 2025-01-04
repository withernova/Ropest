using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public enum InteractiveType
{
    Contact, //接触持续
    Trigger, //手动触发
    Last //停留触发
}

public class InteractiveCaster : MonoBehaviour
{
    public float interactiveRadios;
    public LayerMask interactiveLayerMask;
    public InteractiveType interactiveType;
    public Func<InteractiveBase, bool> condition;
    public float interactiveCD;
    bool inCD;
    public float lastTriggerTime;
    Timer tickCD;

    private void Start()
    {
        tickCD = TimerManager.instance.CreateAndStartTimer(interactiveCD, 1, () => inCD = false);
    }

    private void Update()
    {
        //var a = Physics.OverlapSphere(transform.position, interactiveRadios, interactiveLayerMask);
        if(interactiveType == InteractiveType.Contact)
        {
            ConstantInteractive(Time.deltaTime);
        }

        if(interactiveType == InteractiveType.Last)
        {
            LastInteractive();
        }

        //高亮物体
        //UIManager.instance.showInteractiveUI(a.Length > 0);
    }

    public void TriggerInteractive()
    {
        if (inCD)return;
        
        var a = Physics.OverlapSphere(transform.position, interactiveRadios, interactiveLayerMask);

        foreach (var iTarget in a)
        {
            var interObj = iTarget.GetComponent<InteractiveBase>();
            if (null != interObj)
            {
                if (null != condition && !condition.Invoke(interObj)) continue;
                iTarget.GetComponent<InteractiveBase>().OnStart();
            }
        }

        inCD = true;
        tickCD.Start();
    }

    public void ConstantInteractive(float deltaTime)
    {
        var a = Physics.OverlapSphere(transform.position, interactiveRadios, interactiveLayerMask);

        foreach (var iTarget in a)
        {
            var interObj = iTarget.GetComponent<InteractiveBase>();
            if (null != interObj)
            {
                if (null != condition && !condition.Invoke(interObj)) continue;
                interObj.OnUpdate(deltaTime);
            }
        }
    }

    public void LastInteractive()
    {
        if (inCD) return;

        var a = Physics.OverlapSphere(transform.position, interactiveRadios, interactiveLayerMask);

        foreach (var iTarget in a)
        {
            var interObj = iTarget.GetComponent<InteractiveBase>();
            if (null != interObj)
            {
                if (null != condition && !condition.Invoke(interObj)) continue;
                interObj.lastingTime += Time.deltaTime;
                if (interObj.lastingTime > lastTriggerTime)
                {
                    interObj.OnStart();
                    interObj.lastingTime = 0;
                    inCD = true;
                    tickCD.Start();
                }
            }
        }
    }

    private void OnDrawGizmos()
    {
        Gizmos.color = new Color(1, 0, 0, 0.2f);
        Gizmos.DrawSphere(transform.position,interactiveRadios);
    }

}
