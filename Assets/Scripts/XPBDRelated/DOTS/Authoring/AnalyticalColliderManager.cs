using System.Collections.Generic;
using Unity.Collections;
using Unity.Mathematics;
using UnityEngine;

public class AnalyticalColliderManager : MonoBehaviour
{
    public static AnalyticalColliderManager Instance { get; private set; }

    private static readonly List<AnalyticalColliderSource> _sources = new List<AnalyticalColliderSource>();
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
        UpdateCachedData();
    }

    void OnDestroy()
    {
        if (Instance == this)
            Instance = null;
    }

    public static void Register(AnalyticalColliderSource source)
    {
        if (!_sources.Contains(source))
        {
            _sources.Add(source);
            _dirty = true;
        }
    }

    public static void Unregister(AnalyticalColliderSource source)
    {
        _sources.Remove(source);
        _dirty = true;
    }

    static void UpdateCachedData()
    {
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

    public static NativeArray<AnalyticalColliderData> GetColliderDataForJobs(Allocator allocator)
    {
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
