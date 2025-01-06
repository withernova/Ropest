using DG.Tweening;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;

public class MapManager : MonoBehaviour
{
    public static MapManager Instance;
    public List<string> maps;
    public List<LevelMap> previousMap;
    public List<Material> materials;

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
            "Levels/ClambChanllenge",
            "Levels/CubeChanllenge",
            "Levels/LongRunChallenge",
            "Levels/LongRunChallenge2",
            "Levels/Throw_Hard_1_Chanllenge",
            "Levels/ThrowChanllenge",
            "Levels/RotateClambChanllenge",
            "Levels/RotateClambChanllenge2",
            "Levels/RotateClambChanllenge3",
            "Levels/KiteChanllenge",
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
        GameObject map = Instantiate(Resources.Load<GameObject>("Prefabs/" + path));
      
        LevelMap levelMap = map.GetComponent<LevelMap>();
        //TODO: 详细时间实现
        if (previousMap.Count > 0)
        {
            levelMap.SetBeginningPos(previousMap.Last().GetLasPos().x);
        }
        else
        {
            levelMap.SetBeginningPos(1);
        }
        levelMap.Init(id, 100);


        previousMap.Add(levelMap);
    }

    public void DropMap(LevelMap map)
    {
        map.droped = true;
        previousMap.Remove(map);
        while (previousMap.Count() < 3)
        {
            SpawnMap();
        }
        map.transform.DOMoveY(-10, 3f).onComplete += () =>
        {
            Destroy(map.gameObject);
        };
    }
}