using DG.Tweening;
using UnityEngine;

public class MoveTrap : MonoBehaviour
{
    public float distance;
    public float duration;
    public Vector3 direction;
    public float cd;
    public float delay;
    public Ease ease;

    // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Start()
    {
        Sequence seq = DOTween.Sequence();//
        seq.AppendInterval(delay);
        seq.Append(transform.DOMove(transform.position + distance * direction, duration));
        seq.AppendInterval(cd);
        seq.Append(transform.DOMove(transform.position, duration));
        seq.SetLoops(-1).SetEase(ease);
    }
}
