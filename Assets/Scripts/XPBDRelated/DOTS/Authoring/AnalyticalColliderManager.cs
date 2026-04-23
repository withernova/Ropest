using System.Collections.Generic;
using Unity.Collections;
using Unity.Mathematics;
using UnityEngine;

/// <summary>
/// 解析碰撞体管理器 - 收集场景中的碰撞体数据，供XPBD Burst Job使用
/// 挂载到场景中任意GameObject上，自动收集所有注册的碰撞体
/// 
/// 使用方式：
/// 1. 场景中放置此Manager（单例）
/// 2. 碰撞体GameObject挂载 AnalyticalColliderTag 组件自动注册
/// 3. 或者手动调用 Register/Unregister
/// </summary>
public class AnalyticalColliderManager : MonoBehaviour
{
    public static AnalyticalColliderManager Instance { get; private set; }

    // 已注册的碰撞体列表
    private static readonly List<AnalyticalColliderSource> _sources = new List<AnalyticalColliderSource>();

    // 缓存的碰撞体数据（每帧更新）
    private static AnalyticalColliderData[] _cachedData = new AnalyticalColliderData[0];
    private static bool _dirty = true;

    void Awake()
    {
        if (Instance != null && Instance != this)
        {
            Destroy(gameObject);
            return;
        }
        Instance = this;
    }

    void LateUpdate()
    {
        // 每帧更新碰撞体的Transform数据
        UpdateCachedData();
    }

    void OnDestroy()
    {
        if (Instance == this)
            Instance = null;
    }

    /// <summary>
    /// 注册碰撞体源
    /// </summary>
    public static void Register(AnalyticalColliderSource source)
    {
        if (!_sources.Contains(source))
        {
            _sources.Add(source);
            _dirty = true;
        }
    }

    /// <summary>
    /// 注销碰撞体源
    /// </summary>
    public static void Unregister(AnalyticalColliderSource source)
    {
        _sources.Remove(source);
        _dirty = true;
    }

    /// <summary>
    /// 更新缓存数据（每帧调用，更新Transform）
    /// </summary>
    static void UpdateCachedData()
    {
        // 清理已销毁的源
        _sources.RemoveAll(s => s == null);

        if (_cachedData.Length != _sources.Count)
        {
            _cachedData = new AnalyticalColliderData[_sources.Count];
        }

        for (int i = 0; i < _sources.Count; i++)
        {
            _cachedData[i] = _sources[i].GetColliderData();
        }

        _dirty = false;
    }

    /// <summary>
    /// 获取碰撞体数据的NativeArray副本（供Job使用，调用者负责Dispose）
    /// </summary>
    public static NativeArray<AnalyticalColliderData> GetColliderDataForJobs(Allocator allocator)
    {
        // 确保数据是最新的
        if (_dirty || _cachedData.Length != _sources.Count)
        {
            UpdateCachedData();
        }

        if (_cachedData.Length == 0)
        {
            return new NativeArray<AnalyticalColliderData>(0, allocator);
        }

        var result = new NativeArray<AnalyticalColliderData>(_cachedData.Length, allocator);
        result.CopyFrom(_cachedData);
        return result;
    }
}