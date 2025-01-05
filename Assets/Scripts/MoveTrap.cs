using DG.Tweening;
using UnityEngine;

public class MoveTrap : MonoBehaviour
{
    public float distance;
    public float duration;
    public float cd;
    public float lastTime;
    public Ease ease;

    // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Start()
    {
        Sequence seq = DOTween.Sequence();//
        seq.Append(transform.DOMoveX(transform.position.x + distance, duration));
        seq.AppendInterval(cd);
        seq.Append(transform.DOMoveX(transform.position.x, duration));
        seq.AppendInterval(lastTime);
        seq.SetLoops(-1).SetEase(ease);
    }
}
