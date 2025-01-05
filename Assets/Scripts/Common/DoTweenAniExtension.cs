using DG.Tweening;
using System.Collections.Generic;
using Unity.VisualScripting;
using UnityEngine;
public enum AnimationType
{
    None,
    Scale,
    XScale,
    YScale,
    PumpOnce,
}

public static class DoTweenAniExtension
{
    public static Tween DoPageAnimation(this Transform transform, AnimationType aniType, bool open, float value = -1)
    {
        Tween tween = null;
        transform.gameObject.SetActive(true);
        if (aniType == AnimationType.PumpOnce && open)
        {
            transform.localScale = Vector3.zero;
            transform.DOScale(1, 0.3f).onComplete += () => {
                transform.DOScale(0.95f, 0.05f).onComplete += () => {
                    tween = transform.DOScale(1, 0.1f);
                };
            };
        }
        if (aniType == AnimationType.Scale && open)
        {
            transform.localScale = Vector3.zero;
            tween = transform.DOScale(value == -1 ? 1 : value, 0.3f);
        }

        if (aniType == AnimationType.Scale && !open)
        {
            tween = transform.DOScale(value == -1 ? 0 : value, 0.3f);
        }

        if (aniType == AnimationType.XScale && open)
        {
            transform.localScale = new Vector3(0, 1, 0);
            tween = transform.DOScaleX(1, 0.2f);
        }

        if (aniType == AnimationType.XScale && !open)
        {
            tween = transform.DOScaleX(value == -1 ? 0 : value, 0.2f);
        }

        if (aniType == AnimationType.YScale && open)
        {
            transform.localScale = new Vector3(1, 0, 0);
            tween = transform.DOScaleY(value == -1 ? 1 : value, 0.2f);
        }

        if (aniType == AnimationType.YScale && !open)
        {
            tween = transform.DOScaleY(value == -1 ? 0 : value, 0.2f);
        }

        return tween;
    }
}