using System;
using System.Collections;
using System.Collections.Generic;
using Unity.VisualScripting;
using UnityEngine;

public class ResourcesPool : SingletonForMonoBehaviour<ResourcesPool>
{
    private Dictionary<string, ObjectPool> _objectPools;

    private void Awake()
    {
        _objectPools = new Dictionary<string, ObjectPool>();
    }

    public GameObject Load(string url, int preGenerateCount = 0, Transform PoolPraent = null)
    {
        if (_objectPools.ContainsKey(url))
        {
            return _objectPools[url].GetOne();
        }
        else
        {
            NewPool(url, preGenerateCount, PoolPraent);
            return _objectPools[url].GetOne();
        }
    }

    private void NewPool(string url, int preGenerateCount, Transform poolParent)
    {
        var x = new GameObject();
        x.transform.SetParent(poolParent == null ? this.transform : poolParent);
        x.name = url;
        var c = x.GetOrAddComponent<ObjectPool>();
        var template = Resources.Load<GameObject>("Prefabs/" + url);
        if (null == template)
        {
            Debug.LogError("构建池失败，目标路径" + url + "无实际对象");
        }

        c.Init(template, preGenerateCount);
        _objectPools.Add(url, c);
    }

    public void ReturnOne(string url, GameObject obj)
    {
        var p = _objectPools[url];
        if (null != p)
        {
            p.ReturnOne(obj);
        }
    }
}