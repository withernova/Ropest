using Unity.Collections;
using Unity.Entities;
using Unity.Mathematics;
using UnityEngine;

/// <summary>
/// XPBD碰撞检测工具类 - 提供静态方法供碰撞系统调用
/// 使用传统UnityEngine.Physics API进行碰撞查询（兼容场景中的传统Collider）
/// 仅处理与场景物体的碰撞，自碰撞由Burst Job在SubStep内处理
/// </summary>
public static class XPBDCollisionHelper
{
    // 碰撞辅助GameObject（懒初始化）
    private static GameObject _collisionHelper;
    private static CapsuleCollider _capsuleCollider;
    private static SphereCollider _sphereCollider;

    static void EnsureInitialized()
    {
        if (_collisionHelper != null) return;

        _collisionHelper = new GameObject("XPBD_CollisionHelper");
        _collisionHelper.hideFlags = HideFlags.HideAndDontSave;
        _collisionHelper.layer = LayerMask.NameToLayer("Ignore Raycast");

        _capsuleCollider = _collisionHelper.AddComponent<CapsuleCollider>();
        _capsuleCollider.direction = 2; // Z轴方向
        _capsuleCollider.isTrigger = true;

        _sphereCollider = _collisionHelper.AddComponent<SphereCollider>();
        _sphereCollider.isTrigger = true;
        _sphereCollider.enabled = false;
    }

    /// <summary>
    /// Rope场景碰撞检测 - 对每对相邻粒子使用CapsuleCollider + Physics.ComputePenetration
    /// </summary>
    public static void HandleRopeCollision(NativeArray<float3> posArr, NativeArray<float3> prevArr,
        int numPoints, float radius)
    {
        EnsureInitialized();

        _capsuleCollider.enabled = true;
        _sphereCollider.enabled = false;
        _capsuleCollider.radius = radius;

        for (int i = 0; i < numPoints - 1; i++)
        {
            Vector3 p0 = posArr[i];
            Vector3 p1 = posArr[i + 1];

            Vector3 delta = p1 - p0;
            _collisionHelper.transform.position = p0 + delta / 2f;
            _collisionHelper.transform.LookAt(p0);
            _capsuleCollider.height = delta.magnitude + 2f * radius;

            Collider[] overlaps = Physics.OverlapCapsule(p0, p1, radius + 0.02f, ~LayerMask.GetMask("Ignore Raycast"));

            foreach (var collider in overlaps)
            {
                if (collider.isTrigger) continue;
                if (collider.gameObject == _collisionHelper) continue;

                bool overlapped = Physics.ComputePenetration(
                    _capsuleCollider, _collisionHelper.transform.position, _collisionHelper.transform.rotation,
                    collider, collider.transform.position, collider.transform.rotation,
                    out Vector3 direction, out float distance
                );

                if (overlapped)
                {
                    float3 correction = (float3)(distance * direction);
                    posArr[i] += correction;
                    posArr[i + 1] += correction;

                    ApplyFriction(ref posArr, ref prevArr, i, distance, (float3)direction);
                    ApplyFriction(ref posArr, ref prevArr, i + 1, distance, (float3)direction);
                }
            }
        }
    }

    /// <summary>
    /// Cloth场景碰撞检测 - 对每个粒子使用SphereCollider + Physics.ComputePenetration
    /// 与老方案ClothXPBDSolver.PostSolve保持一致的碰撞策略
    /// </summary>
    public static void HandleClothCollision(NativeArray<float3> posArr, NativeArray<float3> prevArr,
        NativeArray<float3> velArr, NativeArray<float> invMassArr, int numParticles, float collisionRadius, float friction, float dt)
    {
        EnsureInitialized();

        _capsuleCollider.enabled = false;
        _sphereCollider.enabled = true;
        _sphereCollider.radius = collisionRadius;

        float oneOverDt = dt > 0f ? 1f / dt : 0f;

        for (int i = 0; i < numParticles; i++)
        {
            if (invMassArr[i] == 0f) continue;

            // 用当前位置做碰撞检测（与老方案一致）
            _collisionHelper.transform.position = (Vector3)posArr[i];

            Collider[] overlaps = Physics.OverlapSphere((Vector3)posArr[i], collisionRadius + 0.02f, ~LayerMask.GetMask("Ignore Raycast"));

            foreach (var collider in overlaps)
            {
                if (collider.isTrigger) continue;
                if (collider.gameObject == _collisionHelper) continue;

                // 更新helper位置到最新（多碰撞体时位置可能已变）
                _collisionHelper.transform.position = (Vector3)posArr[i];

                bool overlapped = Physics.ComputePenetration(
                    _sphereCollider, _collisionHelper.transform.position, _collisionHelper.transform.rotation,
                    collider, collider.transform.position, collider.transform.rotation,
                    out Vector3 direction, out float distance
                );

                if (overlapped)
                {
                    float3 dir = (float3)direction;
                    // 位置修正：推出碰撞体（与老方案完全一致）
                    posArr[i] += distance * dir;

                    // 摩擦处理：与老方案FrictionForNative一致
                    // 计算位移向量
                    float3 disp = posArr[i] - prevArr[i];
                    // 投影到碰撞面（去除法线分量）
                    float3 tangent = disp - math.dot(disp, dir) * dir;
                    float tangentLen = math.length(tangent);

                    if (tangentLen < 0.01f)
                    {
                        // 几乎没有切向运动，完全停止
                        prevArr[i] = posArr[i];
                    }
                    else
                    {
                        // 沿运动方向施加摩擦（与老方案一致：prevPos += normalized(disp) * 0.25 * distance）
                        float dispLen = math.length(disp);
                        if (dispLen > 1e-8f)
                        {
                            prevArr[i] += (disp / dispLen) * friction * distance;
                        }
                    }

                    // 刚体交互
                    if (dt > 0f && collider.TryGetComponent<Rigidbody>(out var rb))
                    {
                        float3 v1 = (posArr[i] - prevArr[i]) * oneOverDt;
                        float3 v2 = (float3)rb.linearVelocity;
                        float m1 = 1f / invMassArr[i];
                        float m2 = rb.mass;
                        float3 v1_ = (v1 * (m1 - m2) + 2 * m2 * v2) / (m1 + m2);
                        float3 v2_ = (v2 * (m2 - m1) + 2 * m1 * v1) / (m1 + m2);

                        prevArr[i] = posArr[i] - v1_ * dt;
                        rb.linearVelocity = v2_;
                    }
                }
            }

            // 碰撞修正后重新计算速度（老方案也是先碰撞再算速度）
            velArr[i] = (posArr[i] - prevArr[i]) * oneOverDt;
        }
    }

    static void ApplyFriction(ref NativeArray<float3> pos, ref NativeArray<float3> prevPos,
        int i, float distance, float3 direction)
    {
        float3 p1 = pos[i];
        float3 p1_ = prevPos[i];

        float3 p1_p1 = p1 - p1_;
        float3 t1 = p1_p1 - math.dot(p1_p1, direction) * direction;

        if (math.length(t1) < 0.01f)
        {
            prevPos[i] = p1;
            return;
        }

        prevPos[i] += math.normalize(p1_p1) * 0.25f * distance;
    }

    public static void Cleanup()
    {
        if (_collisionHelper != null)
        {
            Object.Destroy(_collisionHelper);
            _collisionHelper = null;
        }
    }
}

/// <summary>
/// XPBD碰撞检测系统 - PostSolve后做一次场景碰撞修正
/// 自碰撞由SubStep内的Burst Job处理，这里只处理与场景物体的碰撞
/// </summary>
[UpdateInGroup(typeof(FixedStepSimulationSystemGroup))]
[UpdateAfter(typeof(RopeSimulationSystem))]
[UpdateAfter(typeof(ClothSimulationSystem))]
public partial class XPBDCollisionSystem : SystemBase
{
    private EntityQuery _ropeQuery;
    private EntityQuery _clothQuery;

    protected override void OnCreate()
    {
        _ropeQuery = GetEntityQuery(
            ComponentType.ReadOnly<RopeTag>(),
            ComponentType.ReadOnly<RopeSolverConfig>(),
            ComponentType.ReadWrite<ParticlePosition>(),
            ComponentType.ReadWrite<ParticlePrevPosition>()
        );
        _clothQuery = GetEntityQuery(
            ComponentType.ReadOnly<ClothTag>(),
            ComponentType.ReadOnly<ClothSolverConfig>(),
            ComponentType.ReadWrite<ParticlePosition>(),
            ComponentType.ReadWrite<ParticlePrevPosition>(),
            ComponentType.ReadWrite<ParticleVelocity>(),
            ComponentType.ReadOnly<ParticleInvMass>()
        );
    }

    protected override void OnStopRunning()
    {
        XPBDCollisionHelper.Cleanup();
    }

    protected override void OnUpdate()
    {
        Dependency.Complete();

        // Rope场景碰撞
        {
            var entities = _ropeQuery.ToEntityArray(Allocator.Temp);
            for (int e = 0; e < entities.Length; e++)
            {
                var entity = entities[e];
                var cfg = EntityManager.GetComponentData<RopeSolverConfig>(entity);
                var positions = EntityManager.GetBuffer<ParticlePosition>(entity);
                var prevPositions = EntityManager.GetBuffer<ParticlePrevPosition>(entity);

                if (positions.Length == 0 || cfg.NumPoints <= 0) continue;

                var posArr = positions.Reinterpret<float3>().AsNativeArray();
                var prevArr = prevPositions.Reinterpret<float3>().AsNativeArray();

                XPBDCollisionHelper.HandleRopeCollision(posArr, prevArr, cfg.NumPoints, cfg.Radius);
            }
            entities.Dispose();
        }

        // Cloth场景碰撞
        {
            var entities = _clothQuery.ToEntityArray(Allocator.Temp);
            for (int e = 0; e < entities.Length; e++)
            {
                var entity = entities[e];
                var cfg = EntityManager.GetComponentData<ClothSolverConfig>(entity);
                var positions = EntityManager.GetBuffer<ParticlePosition>(entity);
                var prevPositions = EntityManager.GetBuffer<ParticlePrevPosition>(entity);
                var invMasses = EntityManager.GetBuffer<ParticleInvMass>(entity);
                var velocities = EntityManager.GetBuffer<ParticleVelocity>(entity);

                if (positions.Length == 0 || cfg.NumParticles <= 0) continue;

                var posArr = positions.Reinterpret<float3>().AsNativeArray();
                var prevArr = prevPositions.Reinterpret<float3>().AsNativeArray();
                var invMassArr = invMasses.Reinterpret<float>().AsNativeArray();
                var velArr = velocities.Reinterpret<float3>().AsNativeArray();

                // 碰撞修正（与老方案一致：一次碰撞修正 + 摩擦 + 速度更新）
                XPBDCollisionHelper.HandleClothCollision(posArr, prevArr, velArr, invMassArr,
                    cfg.NumParticles, cfg.CollisionRadius, cfg.Friction, SystemAPI.Time.DeltaTime);
            }
            entities.Dispose();
        }
    }
}
