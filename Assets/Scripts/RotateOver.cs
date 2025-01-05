using DG.Tweening;
using System.Collections.Generic;
using UnityEngine;

public class RotateOver : MonoBehaviour
{
    public List<Vector3> angles;
    public List<float> cds;
    public List<float> times;

    // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Start()
    {
        Sequence sequence = DOTween.Sequence();
        for (int i = 0; i < angles.Count; i++)
        {
            sequence.Append(transform.DOShakeRotation(1, 5, 5));
            sequence.Append(transform.DOLocalRotate(angles[i], times[i], RotateMode.LocalAxisAdd));
            sequence.AppendInterval(cds[i]);
        }
        sequence.SetLoops(-1);
    }
}
