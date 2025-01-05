using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;
using static UnityEditor.Progress;

public enum InteractiveType
{
    Contact, //接触持续
    Trigger, //手动触发
    Last, //停留触发
    UnDeploy //不负责触发逻辑
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

        var inCondition = a.ToList().Where(item => null != item.GetComponent<InteractiveBase>() && !(null != condition && !condition.Invoke(item.GetComponent<InteractiveBase>())));
        var res = inCondition.Where(item => (transform.position - item.transform.position).magnitude == a.Min(item => (transform.position - item.transform.position).magnitude));
        

        res.First().GetComponent<InteractiveBase>().OnStart();

        inCD = true;
        tickCD.Start();
    }

    public T TriggerInteractiveUndeploy<T>() where T : InteractiveBase
    {
        if (inCD) return null;

        var a = Physics.OverlapSphere(transform.position, interactiveRadios, interactiveLayerMask, QueryTriggerInteraction.Collide);

        var inCondition = a.ToList().Where(item => null != item.GetComponent<T>() && !(null != condition && !condition.Invoke(item.GetComponent<T>())));
        var res = inCondition.ToList().OrderBy(item => (item.transform.position - transform.position).magnitude);
        //Debug.Log(res.Count());
        //Debug.Log(inCondition.Count() + "a" + res.Count());
        if(res.Count() > 0 && res.First().GetComponent<InteractiveBase>().radius != 0 && (transform.position - res.First().transform.position).magnitude < res.First().GetComponent<InteractiveBase>().radius)
        {
            inCD = true;
            tickCD.Start();
            return res.First().GetComponent<T>();
        }

        return null;
    }

    public void ConstantInteractive(float deltaTime)
    {
        var a = Physics.OverlapSphere(transform.position, interactiveRadios, interactiveLayerMask);

        var inCondition = a.ToList().Where(item => null != item.GetComponent<InteractiveBase>() && !(null != condition && !condition.Invoke(item.GetComponent<InteractiveBase>())));
        var res = inCondition.Where(item => (transform.position - item.transform.position).magnitude == a.Min(item => (transform.position - item.transform.position).magnitude));

        res.First().GetComponent<InteractiveBase>().OnUpdate(deltaTime);
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
