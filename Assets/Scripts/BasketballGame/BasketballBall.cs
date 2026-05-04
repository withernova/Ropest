using Unity.Entities;
using Unity.Mathematics;
using UnityEngine;

namespace BasketballGame
{
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

        public void SetInitialVelocity(Vector3 velocity)
        {
            _initialVelocity = velocity;
        }

        void LateUpdate()
        {
            // XPBD 物理暂停时冻结所有游戏逻辑：
            // - 不累加 _lifeTimer，避免暂停中自动销毁；
            // - 不触发生命周期销毁 / killY 检查；
            // - 不更新 Center（物理也没推进，没必要刷新）。
            // 只保留"应用初速度"的兜底——极端情况下玩家一暂停就按 N，
            // 此时球粒子还没应用初速度，下面 TryApplyInitialVelocity 仍需运行一次。
            if (XPBDDebugController.IsPaused)
            {
                if (!_initialVelocityApplied)
                {
                    if (TryApplyInitialVelocity())
                        _initialVelocityApplied = true;
                }
                return;
            }

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
