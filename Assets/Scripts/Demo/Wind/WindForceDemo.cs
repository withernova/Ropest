using Unity.Collections;
using Unity.Entities;
using Unity.Mathematics;
using UnityEngine;

namespace Ropest.Demo
{
    /// <summary>
    /// 录屏专用的"伪风力"驱动器（拖曳模型，稳定不发散）。
    ///
    /// 之前版本是"每 PulseInterval 秒往速度上加一个固定增量"，问题：
    ///   - 不衰减、不封顶 → 几次脉冲后速度爆掉 → XPBD 解算 NaN → Mesh AABB Invalid。
    ///
    /// 现在改为"每帧把粒子速度往目标风速 v_target 软拉"：
    ///     dv = (v_target - v) * (1 - exp(-k * dt))     // 指数趋近
    ///   并对叠加后的速度做模长 clamp，再加一个 yaw sin 摆动让"风"看起来不死板。
    ///
    /// 这样的好处：
    ///   - 粒子速度有理论上限 = |v_target| + |v_perp_clamp|，不会发散；
    ///   - 不需要每帧 Complete + 写 Buffer 这种重操作放在 Update 里循环——同样在 Update 里写，
    ///     但每帧只做一次轻量化加权，不会和 ClothSimulationSystem 抢；
    ///   - 视觉上是"风把布料往一个方向推"的稳态效果，比脉冲更像旗帜。
    /// </summary>
    [DisallowMultipleComponent]
    public class WindForceDemo : MonoBehaviour
    {
        [Tooltip("风的目标速度方向（世界空间，会归一化）")]
        public Vector3 Direction = new Vector3(1f, 0f, 0f);

        [Tooltip("目标风速，米/秒。粒子速度会被软拉向这个目标值。")]
        public float TargetSpeed = 4.0f;

        [Tooltip("拉拽响应速率：值越大风作用越强。建议 1.5~4")]
        [Range(0.1f, 10f)] public float ResponseRate = 2.5f;

        [Tooltip("方向左右摆动幅度（度），0=一直直吹")]
        public float YawWiggleDeg = 25f;

        [Tooltip("摆动频率（Hz）")]
        public float WiggleFrequency = 1.7f;

        [Tooltip("风对高处粒子更敏感的程度。0=均匀，1=只吹顶端")]
        [Range(0f, 1f)] public float HeightBias = 0.4f;

        [Tooltip("粒子速度模长上限，防数值发散。建议 = TargetSpeed * 2.5")]
        public float MaxSpeed = 10f;

        [Tooltip("每帧应用一次（false=按 Interval 间隔应用，省 CPU 但效果略一抖一抖）")]
        public bool ApplyEveryFrame = true;

        [Tooltip("ApplyEveryFrame=false 时的应用间隔（秒）")]
        public float Interval = 0.1f;

        float _timer;
        float _phase;

        void Update()
        {
            float dt = Time.deltaTime;
            _phase += dt;

            if (!ApplyEveryFrame)
            {
                _timer += dt;
                if (_timer < Interval) return;
                dt = _timer;        // 把累积的 dt 一次性补给指数趋近，效果一致
                _timer = 0f;
            }

            ApplyWind(dt);
        }

        void ApplyWind(float dt)
        {
            var world = World.DefaultGameObjectInjectionWorld;
            if (world == null || !world.IsCreated) return;
            var em = world.EntityManager;

            // 在写 ParticleVelocity 前必须先 Complete 所有未完成的物理 Job，
            // 否则会被 DOTS safety check 拦截。
            em.CompleteAllTrackedJobs();

            var query = em.CreateEntityQuery(
                ComponentType.ReadOnly<ClothTag>(),
                ComponentType.ReadWrite<ParticleVelocity>(),
                ComponentType.ReadOnly<ParticlePosition>(),
                ComponentType.ReadOnly<ParticleInvMass>());
            var entities = query.ToEntityArray(Allocator.Temp);

            // 当前帧风向（绕 Y 轴正弦摆动）
            Vector3 dirN = Direction.sqrMagnitude > 1e-6f ? Direction.normalized : Vector3.right;
            float yawDeg = Mathf.Sin(_phase * WiggleFrequency * Mathf.PI * 2f) * YawWiggleDeg;
            Vector3 dir = Quaternion.AngleAxis(yawDeg, Vector3.up) * dirN;
            float3 vTarget = (float3)(dir * TargetSpeed);

            // 指数趋近因子：alpha = 1 - exp(-k*dt)，dt 越大或 k 越大，alpha 越接近 1
            float alpha = 1f - math.exp(-math.max(0f, ResponseRate) * math.max(1e-4f, dt));

            for (int e = 0; e < entities.Length; e++)
            {
                var entity = entities[e];
                var velBuf = em.GetBuffer<ParticleVelocity>(entity);
                var posBuf = em.GetBuffer<ParticlePosition>(entity, true);
                var massBuf = em.GetBuffer<ParticleInvMass>(entity, true);

                int n = velBuf.Length;
                if (n == 0) continue;

                // 计算 Y 范围用于高度权重
                float yMin = float.PositiveInfinity, yMax = float.NegativeInfinity;
                for (int i = 0; i < n; i++)
                {
                    float y = posBuf[i].Value.y;
                    if (math.isfinite(y))
                    {
                        if (y < yMin) yMin = y;
                        if (y > yMax) yMax = y;
                    }
                }
                if (!math.isfinite(yMin) || !math.isfinite(yMax))
                {
                    // 已经是 NaN 帧，不再注入能量
                    continue;
                }
                float yRange = math.max(1e-3f, yMax - yMin);

                for (int i = 0; i < n; i++)
                {
                    if (massBuf[i].Value <= 0f) continue; // 固定点不动
                    float yNorm = (posBuf[i].Value.y - yMin) / yRange;        // 0..1
                    float weight = math.lerp(1f - HeightBias, 1f, yNorm);
                    float a = math.clamp(alpha * weight, 0f, 1f);

                    float3 v = velBuf[i].Value;
                    if (!math.isfinite(v.x) || !math.isfinite(v.y) || !math.isfinite(v.z))
                    {
                        // 单粒子已 NaN，重置，避免污染整体。这种情况一般来源于上一帧已经发散，
                        // 正常路径下不会进入。
                        v = float3.zero;
                    }

                    // 软拉向目标速度（仅水平分量；y 用更小的系数，避免风把布料整体向上吹飞）
                    v.x = math.lerp(v.x, vTarget.x, a);
                    v.z = math.lerp(v.z, vTarget.z, a);
                    v.y = math.lerp(v.y, vTarget.y, a * 0.3f);

                    // 速度模长 clamp
                    float speed2 = math.lengthsq(v);
                    if (speed2 > MaxSpeed * MaxSpeed)
                    {
                        v *= MaxSpeed / math.sqrt(math.max(1e-8f, speed2));
                    }

                    velBuf[i] = new ParticleVelocity { Value = v };
                }
            }

            entities.Dispose();
            query.Dispose();
        }
    }
}
