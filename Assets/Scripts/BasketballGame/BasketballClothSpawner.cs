using System.Collections;
using UnityEngine;

namespace BasketballGame
{
    /// <summary>
    /// 布料目标生成器：周期性地在目标区域内的随机位置生成一张布料，
    /// 由 BasketballGameManager 限制同时存在的数量（超过上限会自动销毁最早的）。
    ///
    /// 布料结构：
    ///  - 通过 ClothRuntimeSpawner 创建
    ///  - 默认顶部两角固定（fixedVertices 设置为 [0, subdivision]，上边缘左右两端）
    ///  - meshOrigin 使用"地下起点"，BasketballClothTarget 内部通过移动固定点实现"升起"
    /// </summary>
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
        public float minDistanceFromPlayer = 6f;
        [Tooltip("X方向随机范围（相对 targetAreaCenter 局部）")]
        public float xSpread = 5f;
        [Tooltip("Z方向额外偏移范围")]
        public float zSpread = 3f;

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

        /// <summary>
        /// 手动生成一张布料（也可被外部调用）
        /// </summary>
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

            // 关键：布料的 "meshOrigin" 是布料左下角（生成用），布料中心应当在 targetOrigin
            // 所以先算布料中心偏移：宽度 length 在 X，高度 width 在 Z
            // 但我们希望布料是竖直的（门帘式）：在 Spawner 生成后，布料是躺平在XZ平面的
            // 为了让它"站起来"，我们让 prefab 本身的 ClothRuntimeSpawner 已经配置为竖直放置（通过 length 指 X，width 指 Y 方向）
            // 实际上 ClothRuntimeSpawner 的 CreateClothMesh 把顶点放在 XZ 平面。
            // 为了简单起见，我们让布料"门帘"的法线朝向玩家，X沿玩家水平方向，Z方向=垂直于地面
            // 这需要改动 ClothRuntimeSpawner，或者我们在 BasketballClothTarget 里把粒子旋转90度。
            // 更简单的方案：给 clothPrefab 配置成"躺在地上"的布，垂直升起 clothRiseHeight 后在高空平铺，
            //              球从下方穿过判定 (球心进入一个以布料中心为心，半径=scoreTriggerRadius 的球体)。
            // 这样"球投进布料"= 球从下/侧飞到半空中这张平铺布料的圆心附近。

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
            go.SetActive(true); // 激活后触发 ClothRuntimeSpawner.Start() 生成粒子

            gm.RegisterCloth(target);
            return target;
        }
    }
}
