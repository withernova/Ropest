using System.Collections.Generic;
using Unity.Collections;
using Unity.Mathematics;
using UnityEngine;

public class AnalyticalColliderSource : MonoBehaviour
{
    [Header("碰撞体类型")]
    public AnalyticalColliderType colliderType = AnalyticalColliderType.Sphere;

    [Header("球体参数")]
    public float sphereRadius = 0.5f;

    [Header("Box参数")]
    public Vector3 boxHalfExtents = new Vector3(0.5f, 0.5f, 0.5f);

    [Header("自动检测")]
    public bool autoDetect = true;

    void OnEnable()
    {
        if (autoDetect)
        {
            DetectCollider();
        }
        AnalyticalColliderManager.Register(this);
    }

    void OnDisable()
    {
        AnalyticalColliderManager.Unregister(this);
    }

    void DetectCollider()
    {
        Vector3 lossy = transform.lossyScale;
        if (TryGetComponent<SphereCollider>(out var sphere))
        {
            colliderType = AnalyticalColliderType.Sphere;
            float maxScale = Mathf.Max(Mathf.Abs(lossy.x),
                Mathf.Max(Mathf.Abs(lossy.y), Mathf.Abs(lossy.z)));
            sphereRadius = sphere.radius * maxScale;
        }
        else if (TryGetComponent<BoxCollider>(out var box))
        {
            colliderType = AnalyticalColliderType.Box;
            boxHalfExtents = new Vector3(
                box.size.x * 0.5f * Mathf.Abs(lossy.x),
                box.size.y * 0.5f * Mathf.Abs(lossy.y),
                box.size.z * 0.5f * Mathf.Abs(lossy.z)
            );
        }
    }

    public AnalyticalColliderData GetColliderData()
    {
        Vector3 lossy = transform.lossyScale;
        Quaternion rot = transform.rotation;

        var data = new AnalyticalColliderData
        {
            Type = colliderType,
            Center = transform.position,
            Rotation = rot,
            InvRotation = math.inverse(rot)
        };

        if (colliderType == AnalyticalColliderType.Sphere)
        {
            if (autoDetect && TryGetComponent<SphereCollider>(out var sphere))
            {
                data.Center = (float3)(transform.position + rot * Vector3.Scale(sphere.center, lossy));
                float maxScale = Mathf.Max(Mathf.Abs(lossy.x),
                    Mathf.Max(Mathf.Abs(lossy.y), Mathf.Abs(lossy.z)));
                data.Radius = sphere.radius * maxScale;
            }
            else
            {
                data.Radius = sphereRadius;
            }
        }
        else
        {
            if (autoDetect && TryGetComponent<BoxCollider>(out var box))
            {
                data.Center = (float3)(transform.position + rot * Vector3.Scale(box.center, lossy));
                data.HalfExtents = new float3(
                    box.size.x * 0.5f * Mathf.Abs(lossy.x),
                    box.size.y * 0.5f * Mathf.Abs(lossy.y),
                    box.size.z * 0.5f * Mathf.Abs(lossy.z)
                );
            }
            else
            {
                data.HalfExtents = boxHalfExtents;
            }
        }

        return data;
    }
}
