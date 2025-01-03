using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ObjectPool : MonoBehaviour
{
    private List<GameObject> _GameObjects;
    private Queue<GameObject> _freeObjects;
    private GameObject _template;

    public void Init(GameObject template, int preGenerateCount = 0)
    {
        _GameObjects = new List<GameObject>();
        _freeObjects = new Queue<GameObject>();

        if (null == template)
        {
            Debug.LogError("这里怎么给了个空资源呢？？");
            return;
        }

        _template = template;
        for (int i = 0; i < preGenerateCount; i++)
        {
            NewOne();
        }
    }

    public GameObject GetOne()
    {
        if (_freeObjects.Count <= 0)
        {
            NewOne();
        }

        var x = _freeObjects.Dequeue();
        x.name = gameObject.name;
        //        x.GetComponent<HisaoMono>().OnPoolObjectCreate();
        x.SetActive(true);
        return x;
    }

    public void ReturnOne(GameObject obj)
    {
        _freeObjects.Enqueue(obj);
        var x = _freeObjects.Dequeue();
        // x.GetComponent<HisaoMono>().OnPoolObjectReturn();
        obj.SetActive(false);
    }

    public void NewOne()
    {
        var x = Instantiate(_template, transform);
        x.SetActive(false);
        _GameObjects.Add(x);
        _freeObjects.Enqueue(x);
    }
}