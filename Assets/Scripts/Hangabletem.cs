using NUnit.Framework;
using System.Collections.Generic;
using UnityEngine;
using static UnityEditor.Progress;

public class Hangabletem : MonoBehaviour
{
    public List<Vector3> handlers;

    private void Start()
    {
        foreach (Transform child in transform)
        {
            handlers.Add(child.position);
        }
    }

    public Vector3 GetTargetPos(Vector3 origin)
    {
        float distance = 1e5f;
        Vector3 ret = handlers[0];
        foreach (Vector3 item in handlers)
        {
            var dis = (item - origin).magnitude;
            if (dis < distance)
            {
                distance = dis;
                ret = item;
            }
        }

        return ret;

    }

}
