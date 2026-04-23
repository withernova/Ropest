using Unity.Entities;
using Unity.Mathematics;
using UnityEngine;

namespace BasketballGame
{
    /// <summary>
    /// 投篮用的"球"对象。挂在球 Prefab 上，伴随 SoftBodyRuntimeSpawner 使用。
    ///
    /// 职责：
    /// 1. 等待 SoftBodyRuntimeSpawner 把 Entity 创建完后，设置初始速度
    /// 2. 每帧读取 ECS 侧粒子位置，计算球质心（用作进球判定、越界检测）
    /// 3. 生命周期：掉到地面下 / 静止过久 / 超时 -> 自毁
    /// </summary>
    [RequireComponent(typeof(SoftBodyRuntimeSpawner))]
    public class BasketballBall : MonoBehaviour
    {
        [Header("生命周期")]
        [Tooltip("最大存活时间（秒），避免球永久留存")]
        public float maxLifetime = 12f;
        [Tooltip("掉到此高度以下立即销毁")]
        public float killY = -20f;

        // 运行时
        private SoftBodyRuntimeSpawner _spawner;
        private EntityManager _em;
        private float _lifeTimer;
        private Vector3 _initialVelocity;
        private bool _initialVelocityApplied;

        public Vector3 Center { get; private set; }

        void Awake()
        {
            _spawner = GetComponent<SoftBodyRuntimeSpawner>();
            if (World.DefaultGameObjectInjectionWorld != null)
            {
                _em = World.DefaultGameObjectInjectionWorld.EntityManager;
            }
        }

        /// <summary>
        /// 由 Shooter 在 Instantiate 后立刻调用，稍后在 LateUpdate 中应用
        /// （因为 Entity 是在 SoftBodyRuntimeSpawner.Start 中创建的）
        /// </summary>
        public void SetInitialVelocity(Vector3 velocity)
        {
            _initialVelocity = velocity;
        }

        void LateUpdate()
        {
            _lifeTimer += Time.deltaTime;

            // 应用初速度：等到 SoftBodyRuntimeSpawner 创建完 Entity 后再写入
            if (!_initialVelocityApplied)
            {
                if (TryApplyInitialVelocity())
                    _initialVelocityApplied = true;
            }

            // 更新球质心（质心 = 所有粒子位置平均）
            UpdateCenter();

            // 生命周期检查
            if (_lifeTimer > maxLifetime || Center.y < killY)
            {
                var gm = BasketballGameManager.Instance;
                if (gm != null) gm.UnregisterBall(this);
                Destroy(gameObject);
            }
        }

        private bool TryApplyInitialVelocity()
        {
            var entity = FindEntityBySelf();
            if (entity == Entity.Null) return false;

            if (!_em.HasBuffer<ParticlePosition>(entity)) return false;

            var posBuf = _em.GetBuffer<ParticlePosition>(entity);
            var prevBuf = _em.GetBuffer<ParticlePrevPosition>(entity);
            var velBuf = _em.GetBuffer<ParticleVelocity>(entity);
            var massBuf = _em.GetBuffer<ParticleInvMass>(entity);

            if (posBuf.Length == 0) return false;

            // 1. 先把粒子整体平移到 transform.position（因为 Spawner 是以 meshOrigin 生成顶点，默认位于世界原点附近）
            //    SpawnPos 已经包含在 transform.position 中
            //    我们假设 Spawner 生成的球中心在局部坐标 meshOrigin + size/2 处；这里直接把所有粒子再加 transform.position
            float3 offset = transform.position;
            for (int i = 0; i < posBuf.Length; i++)
            {
                float3 newPos = posBuf[i].Value + offset;
                posBuf[i] = new ParticlePosition { Value = newPos };
                prevBuf[i] = new ParticlePrevPosition { Value = newPos };
            }

            // 2. 给所有粒子设置相同初速度（XPBD 中会在下一个 sub-step 生效）
            float3 v = _initialVelocity;
            for (int i = 0; i < velBuf.Length; i++)
            {
                velBuf[i] = new ParticleVelocity { Value = v };
            }

            // 3. 确保所有粒子非固定（invMass != 0）
            for (int i = 0; i < massBuf.Length; i++)
            {
                if (massBuf[i].Value <= 0f)
                {
                    massBuf[i] = new ParticleInvMass { Value = 1f };
                }
            }

            return true;
        }

        private void UpdateCenter()
        {
            var entity = FindEntityBySelf();
            if (entity == Entity.Null) return;
            if (!_em.HasBuffer<ParticlePosition>(entity)) return;
            var posBuf = _em.GetBuffer<ParticlePosition>(entity);
            if (posBuf.Length == 0) return;

            float3 sum = float3.zero;
            for (int i = 0; i < posBuf.Length; i++)
            {
                sum += posBuf[i].Value;
            }
            Center = (Vector3)(sum / posBuf.Length);
        }

        /// <summary>
        /// 通过 ManagedMeshReference 反查当前 GameObject 所对应的 Entity。
        /// SoftBodyRuntimeSpawner 会把 gameObject.GetComponent<MeshFilter>() 存到 ManagedMeshReference 上。
        /// </summary>
        private Entity _cachedEntity = Entity.Null;
        private Entity FindEntityBySelf()
        {
            if (_cachedEntity != Entity.Null && _em.Exists(_cachedEntity))
                return _cachedEntity;

            var mf = GetComponent<MeshFilter>();
            if (mf == null) return Entity.Null;

            var query = _em.CreateEntityQuery(
                ComponentType.ReadOnly<SoftBodyTag>(),
                ComponentType.ReadOnly<ManagedMeshReference>());
            var entities = query.ToEntityArray(Unity.Collections.Allocator.Temp);
            try
            {
                foreach (var e in entities)
                {
                    var mref = _em.GetComponentData<ManagedMeshReference>(e);
                    if (mref != null && mref.MeshFilter == mf)
                    {
                        _cachedEntity = e;
                        return e;
                    }
                }
            }
            finally { entities.Dispose(); }
            return Entity.Null;
        }
    }
}
