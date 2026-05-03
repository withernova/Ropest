using Unity.Entities;
using Unity.Mathematics;
using UnityEngine;

namespace BasketballGame
{
    [RequireComponent(typeof(ClothRuntimeSpawner))]
    public class BasketballClothTarget : MonoBehaviour
    {
        public enum State { Rising, Active, Scored, Falling }

        [Header("升降动画参数")]
        [Tooltip("升起到目标高度所用时间")]
        public float riseDuration = 1.2f;
        [Tooltip("投中后落下所用时间")]
        public float fallDuration = 1.5f;
        [Tooltip("落下后再额外停留的销毁延迟")]
        public float fallDestroyDelay = 0.5f;
        [Tooltip("最长存在时间（秒），无论是否投中都会消失")]
        public float maxLifetime = 25f;

        [Header("形状")]
        [Tooltip("四角朝中心收近的距离（米），值越大布料越像个漏斗/袋子")]
        public float cornerPinchAmount = 1f;
        [Tooltip("四角是否朝向玩家方向竖直放置（否则保持水平平铺）")]
        public bool verticalFacingPlayer = true;
        [Tooltip("布面朝玩家方向的仰角（度）。>0 时布面形成\"底边靠近玩家、顶边远离玩家\"的后仰斜盆，开口朝玩家")]
        [Range(0f, 80f)]
        public float tiltAngleTowardsPlayer = 50f;
        [Tooltip("是否在 Init 时随机倾斜角度（启用时每张布料角度在 [tiltAngleMin, tiltAngleMax] 间随机）")]
        public bool randomizeTiltAngle = true;
        [Tooltip("随机倾斜角度下限")]
        [Range(0f, 80f)]
        public float tiltAngleMin = 45f;
        [Tooltip("随机倾斜角度上限")]
        [Range(0f, 80f)]
        public float tiltAngleMax = 60f;
        [Tooltip("四角朝玩家方向凸出的距离（米）。让四角往前凸、中心相对下凹，形成真正能兜住球的'兜'")]
        [Range(0f, 3f)]
        public float cornerPushBackAmount = 0f;
        [Tooltip("把四角整体相对中心抬高的距离（米）。与 cornerPushBackAmount 配合，让四角形成兜的边沿")]
        [Range(0f, 3f)]
        public float cornerLiftUpAmount = 0.9f;

        // 运行时
        private ClothRuntimeSpawner _spawner;
        private EntityManager _em;
        private bool _emReady;
        private Entity _cachedEntity = Entity.Null;
        private Vector3 _startWorldOrigin;
        private Vector3 _targetWorldOrigin;
        private Vector3 _facePlayerDir = Vector3.back; // 布面朝向玩家的方向（单位向量）
        private float _animTimer;
        private float _lifeTimer;
        private float _stayTimer;

        private bool _particleInitialized = false;
        private bool _meshRendererReadyToEnable = false; // InitializeParticleLayout 成功后延迟一帧启用 MeshRenderer，避免闪帧
        private Vector3[] _particleTargetPositions; // 所有粒子的最终世界位置

        public State CurrentState { get; private set; } = State.Rising;

        public Vector3 Center { get; private set; }

        public float StayProgress01
        {
            get
            {
                var gm = BasketballGameManager.Instance;
                if (gm == null || gm.requiredStayDuration <= 0f) return 0f;
                return Mathf.Clamp01(_stayTimer / gm.requiredStayDuration);
            }
        }

        public void Init(Vector3 targetCenterWorld, float visualRiseHeight = 6f)
        {
            _targetWorldOrigin = targetCenterWorld;

            // 随机倾斜角度：每张布料开口角度略有不同，让游戏更有变化
            if (randomizeTiltAngle)
            {
                float lo = Mathf.Min(tiltAngleMin, tiltAngleMax);
                float hi = Mathf.Max(tiltAngleMin, tiltAngleMax);
                tiltAngleTowardsPlayer = UnityEngine.Random.Range(lo, hi);
            }

            // 视觉起点：只是"Rising 阶段展示时的位置偏移起点"，粒子不会真的生成到这里，
            //           所以即使放到地下也不会引发地面碰撞问题。
            _startWorldOrigin = targetCenterWorld + Vector3.down * visualRiseHeight;

            _spawner = GetComponent<ClothRuntimeSpawner>();
            Vector3 cornerOffset = new Vector3(_spawner.length * 0.5f, 0f, _spawner.width * 0.5f);
            // 关键：meshOrigin 以"地下起点"为基准，而非目标世界位置。
            // 否则 ClothRuntimeSpawner.Start() 在第一帧会把一张水平大平板直接渲染到目标高度，
            // 导致画面出现"在相机前闪一下一块水平布"的 bug —— 直到 LateUpdate 里
            // InitializeParticleLayout 跑完，粒子才被旋转/下沉到真正的起点。
            // 放到地下就算 ClothRuntimeSpawner 第一帧渲染了水平初始网格也看不见。
            _spawner.meshOrigin = _startWorldOrigin - cornerOffset;

            // 朝向玩家的方向（XZ 平面）
            var gm = BasketballGameManager.Instance;
            if (gm != null && gm.playerStand != null)
            {
                Vector3 d = gm.playerStand.position - _targetWorldOrigin;
                d.y = 0;
                _facePlayerDir = d.sqrMagnitude > 1e-4f ? d.normalized : Vector3.back;
            }

            Debug.Log($"[ClothTarget/{name}] Init: target={_targetWorldOrigin}, start={_startWorldOrigin}, " +
                      $"meshOrigin={_spawner.meshOrigin}, facePlayerDir={_facePlayerDir}, " +
                      $"fixedVertices.Count={_spawner.fixedVertices?.Count}");
        }

        void Awake()
        {
            _spawner = GetComponent<ClothRuntimeSpawner>();
            if (World.DefaultGameObjectInjectionWorld != null)
            {
                _em = World.DefaultGameObjectInjectionWorld.EntityManager;
                _emReady = true;
                // 包含我们在 LateUpdate 要读/写的所有 Buffer 类型。
                // CompleteDependency() 会把当前挂在这些类型上的仰赖的仿真 Job 全等完，
                // 之后主线程 GetBuffer 就不会触发安全校验异常。
                _simDepQuery = _em.CreateEntityQuery(
                    ComponentType.ReadWrite<ParticlePosition>(),
                    ComponentType.ReadWrite<ParticlePrevPosition>(),
                    ComponentType.ReadWrite<ParticleVelocity>(),
                    ComponentType.ReadOnly<ParticleInvMass>()
                );
            }
        }

        void OnDestroy()
        {
            var gm = BasketballGameManager.Instance;
            if (gm != null) gm.UnregisterCloth(this);
        }

        // 用于在 LateUpdate 中同步仿真依赖的 Query：
        // 覆盖我们在主线程要读/写的 Buffer 类型，对它调用 CompleteDependency() 即可
        // 让 Job System 把这些类型上挂着的仿真 Job 全部等完。
        private EntityQuery _simDepQuery;

        void LateUpdate()
        {
            // 关键：访问仿真 Buffer 前，先等挂在这些类型上的仿真 Job 完成。
            // 仿真系统取消循环内 Complete() 后，主线程必须自己同步，
            // 否则 GetBuffer<ParticlePosition> 会抛 "previously scheduled job writes to ..." 。
            if (_emReady) _simDepQuery.CompleteDependency();

            _lifeTimer += Time.deltaTime;

            UpdateAnimation();
            UpdateCenter();
            UpdateScoring();

            // 兜底：超时强制落下
            if (_lifeTimer > maxLifetime && CurrentState == State.Active)
            {
                FallAndDestroy();
            }
        }

        private void UpdateAnimation()
        {
            var entity = FindEntity();
            if (entity == Entity.Null)
            {
                // 每秒最多 log 1 次避免刷屏
                if ((int)(_lifeTimer * 2f) != (int)((_lifeTimer - Time.deltaTime) * 2f))
                    Debug.LogWarning($"[ClothTarget/{name}] FindEntity returned Null. life={_lifeTimer:F2}s");
                return;
            }

            // 第一次看到 entity：把所有粒子一次性旋转到竖直面向玩家，并记录所有粒子的最终目标位置
            if (!_particleInitialized)
            {
                Debug.Log($"[ClothTarget/{name}] First tick with entity. Calling InitializeParticleLayout...");
                if (!InitializeParticleLayout(entity))
                {
                    Debug.LogWarning($"[ClothTarget/{name}] InitializeParticleLayout returned false (buffer empty?)");
                    return;
                }
                _particleInitialized = true;
                _meshRendererReadyToEnable = true; // 下一帧再启用 MeshRenderer，等 MeshUpdateSystem 把粒子位置刷到 mesh 顶点后再显示
                Debug.Log($"[ClothTarget/{name}] InitializeParticleLayout OK. particleCount={_particleTargetPositions?.Length}, state={CurrentState}");
            }
            else if (_meshRendererReadyToEnable)
            {
                // 现在是 InitializeParticleLayout 成功后的"下一帧"。
                // 此时 ECS 的 MeshUpdateSystem 已经至少跑过一次，clothMesh 顶点已经对应起点布兜位置；
                // 再启用 MeshRenderer 就不会有任何闪帧了。
                var mr = GetComponent<MeshRenderer>();
                if (mr != null) mr.enabled = true;
                _meshRendererReadyToEnable = false;
            }

            bool skipCollision = (CurrentState == State.Rising || CurrentState == State.Falling);
            SetClothCollisionEnabled(entity, !skipCollision);

            switch (CurrentState)
            {
                case State.Rising:
                {
                    _animTimer += Time.deltaTime;
                    float t = Mathf.Clamp01(_animTimer / riseDuration);
                    float s = Mathf.SmoothStep(0f, 1f, t);
                    // 只移动固定点：s=0 → 固定点在 _startWorldOrigin 对应位置；s=1 → 固定点在最终目标位置
                    // 中间粒子由物理仿真带动，自然呈现飘动效果
                    WriteFixedVerticesOnly(entity, s);
                    if ((int)(_animTimer * 5f) != (int)((_animTimer - Time.deltaTime) * 5f))
                        Debug.Log($"[ClothTarget/{name}] Rising t={t:F2} s={s:F2} animTimer={_animTimer:F2}");
                    if (t >= 1f)
                    {
                        CurrentState = State.Active;
                        _animTimer = 0f;
                        Debug.Log($"[ClothTarget/{name}] Rising complete -> Active");
                    }
                    break;
                }
                case State.Active:
                {
                    // 物理正常跑，只钉住固定点在最终目标位置
                    WriteFixedVerticesOnly(entity, 1f);
                    break;
                }
                case State.Falling:
                {
                    _animTimer += Time.deltaTime;
                    float t = Mathf.Clamp01(_animTimer / fallDuration);
                    float s = Mathf.SmoothStep(0f, 1f, t);
                    // 固定点从最终目标位置(s=1)下降到起点位置(s=0)，中间粒子物理仿真
                    WriteFixedVerticesOnly(entity, 1f - s);
                    if (t >= 1f)
                    {
                        Invoke(nameof(DoDestroy), fallDestroyDelay);
                        CurrentState = State.Scored;
                    }
                    break;
                }
            }
        }

        private bool InitializeParticleLayout(Entity entity)
        {
            if (!_em.HasBuffer<ParticlePosition>(entity))
            {
                Debug.LogWarning($"[ClothTarget/{name}] InitializeParticleLayout: entity has no ParticlePosition buffer");
                return false;
            }
            var posBuf = _em.GetBuffer<ParticlePosition>(entity);
            var prevBuf = _em.GetBuffer<ParticlePrevPosition>(entity);
            var velBuf = _em.GetBuffer<ParticleVelocity>(entity);
            if (posBuf.Length == 0)
            {
                Debug.LogWarning($"[ClothTarget/{name}] InitializeParticleLayout: posBuf.Length == 0");
                return false;
            }

            Debug.Log($"[ClothTarget/{name}] InitializeParticleLayout: posBuf.Length={posBuf.Length}, " +
                      $"first particle pos = {posBuf[0].Value}, last = {posBuf[posBuf.Length - 1].Value}");

            // 粒子生成在 _targetWorldOrigin 为中心的水平平面上
            Vector3 originalCenter = _targetWorldOrigin;

            // 计算旋转：把原始法线 +Y 旋到 _facePlayerDir，同时让 +Z（width 方向）变为 +Y（竖直向上）
            Quaternion rot = Quaternion.identity;
            if (verticalFacingPlayer)
            {
                Quaternion r1 = Quaternion.FromToRotation(Vector3.up, _facePlayerDir);
                Vector3 newZ = r1 * Vector3.forward;
                Vector3 projNewZ = Vector3.ProjectOnPlane(newZ, _facePlayerDir);
                Vector3 projUp = Vector3.ProjectOnPlane(Vector3.up, _facePlayerDir);
                if (projNewZ.sqrMagnitude > 1e-6f && projUp.sqrMagnitude > 1e-6f)
                {
                    projNewZ.Normalize();
                    projUp.Normalize();
                    float angle = Vector3.SignedAngle(projNewZ, projUp, _facePlayerDir);
                    Quaternion r2 = Quaternion.AngleAxis(angle, _facePlayerDir);
                    rot = r2 * r1;
                }
                else
                {
                    rot = r1;
                }
            }

            Vector3 tiltAxis = Vector3.Cross(Vector3.up, _facePlayerDir).normalized;
            if (tiltAxis.sqrMagnitude < 1e-6f) tiltAxis = Vector3.right;
            if (verticalFacingPlayer && Mathf.Abs(tiltAngleTowardsPlayer) > 1e-3f)
            {
                // 符号：正角度 → 顶边后仰、底边靠近玩家（相对之前符号反向一次，修正开口方向）
                Quaternion tiltRot = Quaternion.AngleAxis(-tiltAngleTowardsPlayer, tiltAxis);
                rot = tiltRot * rot;
            }

            Debug.Log($"[ClothTarget/{name}] rotation euler={rot.eulerAngles}, " +
                      $"testRotate(+Y)={rot * Vector3.up}, testRotate(+Z)={rot * Vector3.forward}, testRotate(+X)={rot * Vector3.right}");

            // 判断哪些索引是"固定点"
            var fixedSet = new System.Collections.Generic.HashSet<int>();
            if (_spawner.fixedVertices != null)
            {
                foreach (var idx in _spawner.fixedVertices) fixedSet.Add(idx);
            }

            // "兜"方向：朝向玩家的水平方向（让四角朝玩家凸出 → 中心相对朝玩家反方向下凹 → 形成朝玩家张口的盆）
            Vector3 pushBackDir = _facePlayerDir;

            // 计算每个粒子的目标世界位置
            _particleTargetPositions = new Vector3[posBuf.Length];
            for (int i = 0; i < posBuf.Length; i++)
            {
                Vector3 p = (Vector3)posBuf[i].Value;
                Vector3 rel = p - originalCenter;
                Vector3 rotated = rot * rel;
                Vector3 finalPos = originalCenter + rotated;

                // 固定点（四角）形状处理：收拢 + 后推 + 抬高 → 形成兜边
                if (fixedSet.Contains(i))
                {
                    if (cornerPinchAmount > 0f)
                    {
                        Vector3 toCenter = _targetWorldOrigin - finalPos;
                        float dist = toCenter.magnitude;
                        if (dist > 1e-4f)
                        {
                            float pinch = Mathf.Min(cornerPinchAmount, dist * 0.5f);
                            finalPos += toCenter / dist * pinch;
                        }
                    }
                    // 把四角朝玩家方向凸出 + 向上抬，让中心相对朝玩家反方向凹陷，形成朝玩家张口的"盆"
                    finalPos += pushBackDir * cornerPushBackAmount;
                    finalPos += Vector3.up * cornerLiftUpAmount;
                }
                _particleTargetPositions[i] = finalPos;
            }

            Vector3 initialOffset = _startWorldOrigin - _targetWorldOrigin;
            for (int i = 0; i < posBuf.Length; i++)
            {
                Vector3 p = _particleTargetPositions[i] + initialOffset;
                posBuf[i] = new ParticlePosition { Value = (float3)p };
                prevBuf[i] = new ParticlePrevPosition { Value = (float3)p };
                velBuf[i] = new ParticleVelocity { Value = float3.zero };
            }

            Debug.Log($"[ClothTarget/{name}] After init: particle[0] written to {posBuf[0].Value}, " +
                      $"particle[last] = {posBuf[posBuf.Length - 1].Value}, initialOffset={initialOffset}");

            return true;
        }

        private void WriteFixedVerticesOnly(Entity entity, float s)
        {
            if (_particleTargetPositions == null) return;
            if (!_em.HasBuffer<ParticlePosition>(entity)) return;
            var posBuf = _em.GetBuffer<ParticlePosition>(entity);
            var prevBuf = _em.GetBuffer<ParticlePrevPosition>(entity);
            var massBuf = _em.GetBuffer<ParticleInvMass>(entity);

            var fixedList = _spawner.fixedVertices;
            if (fixedList == null) return;

            // 计算固定点的升降偏移
            Vector3 offset = Vector3.Lerp(_startWorldOrigin - _targetWorldOrigin, Vector3.zero, s);

            for (int k = 0; k < fixedList.Count; k++)
            {
                int idx = fixedList[k];
                if (idx < 0 || idx >= posBuf.Length) continue;
                if (idx >= _particleTargetPositions.Length) continue;
                if (massBuf[idx].Value != 0f) continue;
                Vector3 p = _particleTargetPositions[idx] + offset;
                posBuf[idx] = new ParticlePosition { Value = (float3)p };
                prevBuf[idx] = new ParticlePrevPosition { Value = (float3)p };
            }
        }

        private void UpdateCenter()
        {
            var entity = FindEntity();
            if (entity == Entity.Null) return;
            if (!_em.HasBuffer<ParticlePosition>(entity)) return;
            var posBuf = _em.GetBuffer<ParticlePosition>(entity);
            if (posBuf.Length == 0) return;

            float3 sum = float3.zero;
            for (int i = 0; i < posBuf.Length; i++) sum += posBuf[i].Value;
            Center = (Vector3)(sum / posBuf.Length);
        }

        private void UpdateScoring()
        {
            if (CurrentState != State.Active) return;
            var gm = BasketballGameManager.Instance;
            if (gm == null) return;

            // 遍历场上的球，看是否有任意球处于触发半径内
            bool anyBallInside = false;
            var balls = FindObjectsOfType<BasketballBall>();
            float r2 = gm.scoreTriggerRadius * gm.scoreTriggerRadius;
            foreach (var b in balls)
            {
                if (b == null) continue;
                if ((b.Center - Center).sqrMagnitude <= r2)
                {
                    anyBallInside = true;
                    break;
                }
            }

            if (anyBallInside)
            {
                _stayTimer += Time.deltaTime;
                if (_stayTimer >= gm.requiredStayDuration)
                {
                    gm.AddScore(1, this);
                    FallAndDestroy();
                }
            }
            else
            {
                _stayTimer = Mathf.Max(0f, _stayTimer - Time.deltaTime * 0.5f);
            }
        }

        public void FallAndDestroy()
        {
            if (CurrentState == State.Falling || CurrentState == State.Scored) return;
            CurrentState = State.Falling;
            _animTimer = 0f;
        }

        private void DoDestroy()
        {
            Destroy(gameObject);
        }

        private bool _lastCollisionEnabled = true;
        private void SetClothCollisionEnabled(Entity entity, bool enabled)
        {
            if (entity == Entity.Null) return;
            if (!_em.HasComponent<ClothSolverConfig>(entity)) return;
            // 只在状态变化时写，避免每帧无谓 ECS 写入
            if (_lastCollisionEnabled == enabled && _cachedEntity == entity) return;

            var cfg = _em.GetComponentData<ClothSolverConfig>(entity);
            cfg.SkipAnalyticalCollision = !enabled;
            cfg.SkipCrossBodyCollision = !enabled;
            _em.SetComponentData(entity, cfg);
            _lastCollisionEnabled = enabled;
        }

        private Entity FindEntity()
        {
            if (_cachedEntity != Entity.Null && _em.Exists(_cachedEntity))
                return _cachedEntity;

            var mf = GetComponent<MeshFilter>();
            if (mf == null)
            {
                if ((int)(_lifeTimer * 2f) != (int)((_lifeTimer - Time.deltaTime) * 2f))
                    Debug.LogWarning($"[ClothTarget/{name}] FindEntity: no MeshFilter on GameObject yet. life={_lifeTimer:F2}s");
                return Entity.Null;
            }

            var query = _em.CreateEntityQuery(ComponentType.ReadOnly<ClothTag>(),
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
                        Debug.Log($"[ClothTarget/{name}] FindEntity: matched entity={e.Index}:{e.Version}. (Total cloth entities={entities.Length})");
                        return e;
                    }
                }
                if ((int)(_lifeTimer * 2f) != (int)((_lifeTimer - Time.deltaTime) * 2f))
                    Debug.LogWarning($"[ClothTarget/{name}] FindEntity: no match. entities={entities.Length}, myMf={mf.GetInstanceID()}");
            }
            finally { entities.Dispose(); }
            return Entity.Null;
        }
    }
}
