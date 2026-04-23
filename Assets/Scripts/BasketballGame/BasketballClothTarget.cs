using Unity.Entities;
using Unity.Mathematics;
using UnityEngine;

namespace BasketballGame
{
    /// <summary>
    /// 布料目标：一张从地下升起的"布门帘"。
    /// 职责：
    /// 1. 生命周期：升起 -> 等待投中 -> 投中后落下销毁
    /// 2. 检测球是否持续停留在"布料中心 + 半径"球体内，满足时长后计分
    ///
    /// 构造由 BasketballClothSpawner 完成（挂好 ClothRuntimeSpawner 并设置好参数）
    /// </summary>
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
        [Tooltip("四角朝玩家方向凸出的距离（米）。让四角往前凸、中心相对下凹，形成真正能兜住球的'兜'")]
        [Range(0f, 3f)]
        public float cornerPushBackAmount = 0f;
        [Tooltip("把四角整体相对中心抬高的距离（米）。与 cornerPushBackAmount 配合，让四角形成兜的边沿")]
        [Range(0f, 3f)]
        public float cornerLiftUpAmount = 0.9f;

        // 运行时
        private ClothRuntimeSpawner _spawner;
        private EntityManager _em;
        private Entity _cachedEntity = Entity.Null;
        private Vector3 _startWorldOrigin;  // 升起动画"视觉起点"：目标位置下方一段距离（仅视觉用，粒子不会真的生成到这里）
        private Vector3 _targetWorldOrigin; // 布料最终中心位置（粒子实际生成于此，避免穿过地面）
        private Vector3 _facePlayerDir = Vector3.back; // 布面朝向玩家的方向（单位向量）
        private float _animTimer;
        private float _lifeTimer;
        private float _stayTimer;

        // 粒子初始化：ClothRuntimeSpawner.Start 跑完后，第一次 LateUpdate 时
        //   - 把所有粒子从"水平平铺"旋转到"竖直面向玩家"
        //   - 记录每个粒子的"最终世界目标位置"（含四角朝中心收近）
        //
        // 升降动画策略（修复：全程保持物理模拟）：
        //   - Rising/Falling：每帧只把"固定点（四角）"平滑位移到 目标位置 + 纵向偏移，
        //     中间粒子由 XPBD 物理正常仿真（重力 + 距离约束 + 碰撞），
        //     因此布料在升起/落下的过程中会自然飘动，不再是僵硬的整体平移。
        //   - Active：只钉固定点在最终目标位置（与之前一致）。
        //
        // 初始化时：把所有粒子整体放到"Rising 阶段 s=0"位置（地下附近，但仍在空中，
        //   不触发地面碰撞），让物理从那里开始跑；第一帧起固定点就开始往上移动，
        //   物理把中间粒子顺势拉起来。
        private bool _particleInitialized = false;
        private Vector3[] _particleTargetPositions; // 所有粒子的最终世界位置

        public State CurrentState { get; private set; } = State.Rising;

        /// <summary>
        /// 布料世界空间中心（粒子质心），每帧更新
        /// </summary>
        public Vector3 Center { get; private set; }

        /// <summary>
        /// 当前停留计时 0..required，用于UI显示
        /// </summary>
        public float StayProgress01
        {
            get
            {
                var gm = BasketballGameManager.Instance;
                if (gm == null || gm.requiredStayDuration <= 0f) return 0f;
                return Mathf.Clamp01(_stayTimer / gm.requiredStayDuration);
            }
        }

        // =========================================================
        // 初始化：由 Spawner 调用，设定升起起点/终点（布料中心的世界坐标）
        // =========================================================
        public void Init(Vector3 targetCenterWorld, float visualRiseHeight = 6f)
        {
            _targetWorldOrigin = targetCenterWorld;
            // 视觉起点：只是"Rising 阶段展示时的位置偏移起点"，粒子不会真的生成到这里，
            //           所以即使放到地下也不会引发地面碰撞问题。
            _startWorldOrigin = targetCenterWorld + Vector3.down * visualRiseHeight;

            _spawner = GetComponent<ClothRuntimeSpawner>();
            // 关键修复：直接把 meshOrigin 设为"目标位置左下角"，粒子一开始就在空中目标位置生成，
            //          不会穿过地面，物理仿真稳定。
            //          "从地下升起"的视觉效果由 LateUpdate 中覆盖写入粒子位置实现。
            Vector3 cornerOffset = new Vector3(_spawner.length * 0.5f, 0f, _spawner.width * 0.5f);
            _spawner.meshOrigin = _targetWorldOrigin - cornerOffset;

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
            }
        }

        void OnDestroy()
        {
            var gm = BasketballGameManager.Instance;
            if (gm != null) gm.UnregisterCloth(this);
        }

        void LateUpdate()
        {
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

        // =========================================================
        // 升降动画（修复：全程保持物理模拟）：
        //   - Rising/Falling 阶段：只把固定点（四个角）在"地下起点 ↔ 最终目标位置"之间平滑位移，
        //     其余粒子由 XPBD 物理正常仿真，从而布料在升起/落下过程中也会自然飘动；
        //   - Active 阶段：固定点钉在最终目标位置。
        // =========================================================
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
                Debug.Log($"[ClothTarget/{name}] InitializeParticleLayout OK. particleCount={_particleTargetPositions?.Length}, state={CurrentState}");
            }

            // 根据当前状态写入碰撞开关：
            //   - Rising/Falling：关闭和场景（地板）+ 跨体（球）的碰撞，
            //     避免四角强制动画把中间粒子瞬间塞进地面/球，引起 mesh 明显穿模/扭曲；
            //   - Active：开启全部碰撞，恢复正常交互。
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

        /// <summary>
        /// 第一次拿到 Entity 时：
        ///  1) 把所有粒子绕"布料中心"旋转，使布面法线指向玩家（竖直朝向玩家）
        ///  2) 对四角（固定点）按 cornerPinchAmount 朝中心收近
        ///  3) 记录每个粒子"最终世界目标位置"到 _particleTargetPositions
        /// </summary>
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

            // 额外倾斜：让布面绕水平右轴向玩家反方向后仰 tiltAngleTowardsPlayer 度
            // 效果：布面变成"底边近玩家、顶边远玩家"的后仰斜面 + 中心相对玩家反方向下凹 = 朝玩家方向张口的斜盆
            //      球从玩家方向飞来、下落时会打在斜面上被导向凹陷中心，不再弹走
            // 旋转轴 = 与水平面平行、与 _facePlayerDir 垂直的轴（即水平"右向量"）
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

            // 立即把粒子写到"Rising 阶段 s=0"位置（视觉起点：目标位置下方 visualRiseHeight 米）
            // 并清零速度 + prevPos = pos
            // 注意：中间粒子一起偏移到起点，这样它们和固定点一开始就在同一高度，
            //       距离约束不会一瞬间把布料拉得剧烈变形；随后固定点升起时，
            //       物理会自然地把整张布从地下拉到空中，过程中有飘动。
            Vector3 initialOffset = _startWorldOrigin - _targetWorldOrigin; // = (0, -visualRiseHeight, 0)
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

        /// <summary>
        /// 只把固定点（invMass=0 的四个角）平滑位移到"目标位置 + 升降偏移"。
        /// 其余粒子完全由 XPBD 物理仿真接管，因此布料在升起/落下过程中也会自然飘动。
        /// </summary>
        /// <param name="s">升降进度：0=在地下起点，1=在空中目标位置</param>
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

        // =========================================================
        // 计分：球质心持续停留在"布料中心+触发半径"内达到 requiredStayDuration
        // =========================================================
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

        // =========================================================
        // 状态切换
        // =========================================================
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

        // =========================================================
        // 碰撞开关：通过写 ClothSolverConfig 的 Skip 标志，
        // 让 ClothSimulationSystem / CrossBodyCollisionSystem 本帧跳过布料碰撞
        // =========================================================
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

        // =========================================================
        // 实体查找（缓存）
        // =========================================================
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
