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

    void Start()
    {
        Sequence seq = DOTween.Sequence();//
        seq.AppendInterval(delay);
        seq.Append(transform.DOLocalMove(transform.localPosition + distance * direction, duration));
        seq.AppendInterval(cd);
        seq.Append(transform.DOLocalMove(transform.localPosition, duration));
        seq.SetLoops(-1).SetEase(ease);
    }
}
