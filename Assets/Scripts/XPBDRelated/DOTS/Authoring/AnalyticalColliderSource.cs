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
    /// 自动检测Unity Collider并读取参数（仅用作 Inspector 缓存 / autoDetect=false 时的 fallback）
    /// 运行时真实参数由 GetColliderData() 每帧重新计算
    /// </summary>
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

    /// <summary>
    /// 获取当前帧的碰撞体数据
    /// 注意：必须把 Unity Collider 的 center（相对 transform 的局部偏移）
    /// 加到 world Center 上，否则当 collider.center != (0,0,0) 时，
    /// 解析碰撞的位置会和视觉 Mesh 错位（粒子会撞在一个看不见的偏移位置）。
    /// </summary>
    public AnalyticalColliderData GetColliderData()
    {
        Vector3 lossy = transform.lossyScale;
        Quaternion rot = transform.rotation;

        var data = new AnalyticalColliderData
        {
            Type = colliderType,
            // Center 先用 transform.position 打底，下面根据具体 Collider 的 center 再加偏移
            Center = transform.position,
            Rotation = rot,
            InvRotation = math.inverse(rot)
        };

        if (colliderType == AnalyticalColliderType.Sphere)
        {
            if (autoDetect && TryGetComponent<SphereCollider>(out var sphere))
            {
                // 把 SphereCollider.center（local space）变换到 world space 并叠加
                Vector3 localCenter = Vector3.Scale(sphere.center, lossy);
                data.Center = (float3)(transform.position + rot * localCenter);

                // Unity SphereCollider 真实半径 = radius * max(|x|,|y|,|z|) lossyScale
                float maxScale = Mathf.Max(Mathf.Abs(lossy.x),
                    Mathf.Max(Mathf.Abs(lossy.y), Mathf.Abs(lossy.z)));
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
                // 把 BoxCollider.center（local space）变换到 world space 并叠加
                Vector3 localCenter = Vector3.Scale(box.center, lossy);
                data.Center = (float3)(transform.position + rot * localCenter);

                // Unity BoxCollider 真实半尺寸 = size * 0.5 * |lossyScale|（逐分量）
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
