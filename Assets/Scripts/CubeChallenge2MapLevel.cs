using System.Collections.Generic;
using Unity.VisualScripting;
using UnityEngine;

public class CubeChallenge2MapLevel : LevelMap
{
    public List<Transform> toReset = new List<Transform>();
    List<Vector3> value = new List<Vector3>();

    private void Awake()
    {
        toReset.ForEach(t => { value.Add(t.localPosition); });
    }

    public override void Init(int i, int lastTime)
    {
        for (int j = 0; j < toReset.Count; j++)
        {
            Rigidbody rigidbody =  toReset[i].GetComponent<Rigidbody>();
            rigidbody.isKinematic = true;
            toReset[j].rotation = Quaternion.identity;
            toReset[j].position = value[j] += transform.position;
            TimerManager.instance.CreateAndStartTimer(1f, 1, () => rigidbody.GetComponent<Rigidbody>().isKinematic = false);
        }
        base.Init(i, lastTime);
    }
}
