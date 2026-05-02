using System.Collections;
using UnityEngine;

namespace BasketballGame
{
    public class BasketballClothSpawner : MonoBehaviour
    {
        [Header("预制体（挂有 ClothRuntimeSpawner + BasketballClothTarget）")]
        public GameObject clothPrefab;

        [Header("自动生成")]
        [Tooltip("是否在游戏开始时自动持续生成")]
        public bool autoSpawn = true;
        [Tooltip("两次生成之间的最小间隔（秒）")]
        public float spawnInterval = 4.0f;
        [Tooltip("初次延迟（秒）")]
        public float initialDelay = 0.5f;

        [Header("位置随机")]
        [Tooltip("最小距离玩家")]
        public float minDistanceFromPlayer = 4f;
        [Tooltip("X方向随机范围（相对 targetAreaCenter 局部）")]
        public float xSpread = 3.5f;
        [Tooltip("Z方向额外偏移范围")]
        public float zSpread = 2f;

        private Coroutine _loop;

        void Start()
        {
            if (autoSpawn)
            {
                _loop = StartCoroutine(SpawnLoop());
            }
        }

        void OnDisable()
        {
            if (_loop != null)
            {
                StopCoroutine(_loop);
                _loop = null;
            }
        }

        private IEnumerator SpawnLoop()
        {
            yield return new WaitForSeconds(initialDelay);
            while (true)
            {
                var gm = BasketballGameManager.Instance;
                if (gm != null && gm.ActiveClothCount < gm.maxActiveCloths)
                {
                    SpawnOne();
                }
                yield return new WaitForSeconds(spawnInterval);
            }
        }

        public BasketballClothTarget SpawnOne()
        {
            var gm = BasketballGameManager.Instance;
            if (gm == null || clothPrefab == null) return null;

            Vector3 center = gm.targetAreaCenter != null
                ? gm.targetAreaCenter.position
                : transform.position;

            // 随机一个位置：x 在 [-xSpread, xSpread]，z 在 [-zSpread, zSpread]，朝玩家的方向稍微偏一点
            float x = Random.Range(-xSpread, xSpread);
            float z = Random.Range(-zSpread, zSpread);

            Vector3 playerPos = gm.playerStand != null ? gm.playerStand.position : Vector3.zero;
            Vector3 targetOrigin = center + new Vector3(x, gm.clothRiseHeight, z);

            // 保证距离玩家至少 min 距离
            Vector3 flat = targetOrigin - playerPos; flat.y = 0;
            if (flat.magnitude < minDistanceFromPlayer)
            {
                flat = flat.normalized * minDistanceFromPlayer;
                targetOrigin = new Vector3(playerPos.x + flat.x, targetOrigin.y, playerPos.z + flat.z);
            }

            // Instantiate 时不激活，以便先配置参数再启动 Spawner
            var go = Instantiate(clothPrefab);
            go.name = $"ClothTarget_{Time.frameCount}";

            // 读 ClothRuntimeSpawner 确认 prefab 合法
            var spawner = go.GetComponent<ClothRuntimeSpawner>();
            if (spawner == null)
            {
                Debug.LogError("[BasketballClothSpawner] clothPrefab 缺少 ClothRuntimeSpawner");
                Destroy(go);
                return null;
            }

            // 初始化 BasketballClothTarget（传入"布料中心"的目标世界位置，Init 内部会算 meshOrigin）
            var target = go.GetComponent<BasketballClothTarget>();
            if (target == null) target = go.AddComponent<BasketballClothTarget>();
            target.Init(targetOrigin);
            go.SetActive(true);

            gm.RegisterCloth(target);
            return target;
        }
    }
}
