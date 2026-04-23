using System.Collections.Generic;
using Unity.Collections;
using Unity.Mathematics;
using UnityEngine;

/// <summary>
/// 解析碰撞体源 - 挂载到带有碰撞体的GameObject上
/// 自动检测SphereCollider或BoxCollider并注册到管理器
/// </summary>
public class AnalyticalColliderSource : MonoBehaviour
{
    [Header("碰撞体类型（自动检测，也可手动指定）")]
    public AnalyticalColliderType colliderType = AnalyticalColliderType.Sphere;

    [Header("球体参数（当类型为Sphere时）")]
    public float sphereRadius = 0.5f;

    [Header("Box参数（当类型为Box时）")]
    public Vector3 boxHalfExtents = new Vector3(0.5f, 0.5f, 0.5f);

    [Header("自动从Unity Collider读取参数")]
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

    /// <summary>
    /// 自动检测Unity Collider并读取参数
    /// </summary>
    void DetectCollider()
    {
        if (TryGetComponent<SphereCollider>(out var sphere))
        {
            colliderType = AnalyticalColliderType.Sphere;
            // 考虑缩放
            float maxScale = Mathf.Max(transform.lossyScale.x,
                Mathf.Max(transform.lossyScale.y, transform.lossyScale.z));
            sphereRadius = sphere.radius * maxScale;
        }
        else if (TryGetComponent<BoxCollider>(out var box))
        {
            colliderType = AnalyticalColliderType.Box;
            Vector3 scale = transform.lossyScale;
            boxHalfExtents = Vector3.Scale(box.size * 0.5f, scale);
        }
    }

    /// <summary>
    /// 获取当前帧的碰撞体数据
    /// </summary>
    public AnalyticalColliderData GetColliderData()
    {
        var data = new AnalyticalColliderData
        {
            Type = colliderType,
            Center = transform.position,
            Rotation = transform.rotation,
            InvRotation = math.inverse(transform.rotation)
        };

        if (colliderType == AnalyticalColliderType.Sphere)
        {
            // 运行时也考虑缩放变化
            if (autoDetect && TryGetComponent<SphereCollider>(out var sphere))
            {
                float maxScale = Mathf.Max(transform.lossyScale.x,
                    Mathf.Max(transform.lossyScale.y, transform.lossyScale.z));
                data.Radius = sphere.radius * maxScale;
            }
            else
            {
                data.Radius = sphereRadius;
            }
        }
        else // Box
        {
            if (autoDetect && TryGetComponent<BoxCollider>(out var box))
            {
                Vector3 scale = transform.lossyScale;
                data.HalfExtents = Vector3.Scale(box.size * 0.5f, scale);
            }
            else
            {
                data.HalfExtents = boxHalfExtents;
            }
        }

        return data;
    }
}
