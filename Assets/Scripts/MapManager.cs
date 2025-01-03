using DG.Tweening;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;

public class MapManager : MonoBehaviour
{
    public static MapManager Instance;
    public List<string> maps;
    public List<LevelMap> previousMap;

    public void Awake()
    {
        if (Instance == null)
        {
            Instance = this;
        }
    }
    public void Start()
    {
        previousMap = new List<LevelMap>();
        maps = new List<string>()
        {
            "Levels/BarChanllenge",
            "Levels/RowChanllenge",
            "Levels/Level3",
        };

        for (int i = 0; i < 3; i++)
        {
            SpawnMap();
        }
    }

    public void SpawnMap()
    {
        int id = 0;
        do
        {
            id = Random.Range(0, maps.Count);
        } while (previousMap.Any(map => map.id == id));

        string path = maps[id];
        GameObject map = ResourcesPool.Instance.Load(path, 10);
      
        LevelMap levelMap = map.GetComponent<LevelMap>();
        //TODO: 详细时间实现
        levelMap.Init(id, Random.Range(1,10));

        if(previousMap.Count > 0) { 
            levelMap.SetBeginningPos(previousMap.Last().GetLasPos().x);
        }
        previousMap.Add(levelMap);
    }

    public void DropMap(LevelMap map)
    {
        if(map != previousMap.First())
        {
            DropMap(previousMap.First());
        }
        map.droped = true;
        previousMap.Remove(map);
        while (previousMap.Count() < 3)
        {
            SpawnMap();
        }
        map.transform.DOMoveY(-10, 3f).onComplete += () =>
        {
            ResourcesPool.Instance.ReturnOne(maps[map.id], map.gameObject);
        };
    }
}