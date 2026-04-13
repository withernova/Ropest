using Unity.Collections;
using Unity.Entities;
using Unity.Jobs;
using Unity.Mathematics;
using UnityEngine;

/// <summary>
/// Rope Mesh更新系统 - 在Presentation阶段运行
/// 将ECS中的粒子数据回写到Unity Mesh
/// </summary>
[UpdateInGroup(typeof(PresentationSystemGroup))]
public partial class RopeMeshUpdateSystem : SystemBase
{
    private EntityQuery _ropeQuery;

    protected override void OnCreate()
    {
        _ropeQuery = GetEntityQuery(
            ComponentType.ReadOnly<RopeTag>(),
            ComponentType.ReadOnly<RopeSolverConfig>(),
            ComponentType.ReadOnly<ParticlePosition>(),
            ComponentType.ReadOnly<GhostPosition>()
        );
        RequireForUpdate(_ropeQuery);
    }

    protected override void OnUpdate()
    {
        var entities = _ropeQuery.ToEntityArray(Allocator.Temp);

        for (int e = 0; e < entities.Length; e++)
        {
            var entity = entities[e];

            // 获取托管组件
            if (!EntityManager.HasComponent<ManagedMeshReference>(entity)) continue;
            var meshRef = EntityManager.GetComponentObject<ManagedMeshReference>(entity);
            if (meshRef == null || meshRef.Mesh == null) continue;

            var cfg = EntityManager.GetComponentData<RopeSolverConfig>(entity);
            var positions = EntityManager.GetBuffer<ParticlePosition>(entity);
            var ghostPositions = EntityManager.GetBuffer<GhostPosition>(entity);

            // 数据有效性检查
            if (positions.Length == 0 || ghostPositions.Length == 0 || cfg.NumPoints <= 0) continue;

            var posArr = positions.Reinterpret<float3>().AsNativeArray();
            var ghostArr = ghostPositions.Reinterpret<float3>().AsNativeArray();

            int numPoints = cfg.NumPoints;
            int subdivision = cfg.Subdivision;
            float radius = cfg.Radius;

            // 计算输出顶点数量
            int interpolations = 10;
            int numOutputVerts = 1 + (numPoints - 2) * (interpolations - 1) * subdivision + 1;

            // 使用Burst Job进行Cosserat插值
            var outputVerts = new NativeArray<float3>(numOutputVerts, Allocator.TempJob);

            var renderJob = new RopeRenderingJob
            {
                PointPos = posArr,
                GhostPos = ghostArr,
                OutputVertices = outputVerts,
                Subdivision = subdivision,
                Radius = radius,
                NumPoints = numPoints
            };
            renderJob.Schedule(default(JobHandle)).Complete();

            // 回写到Mesh
            var vertices = new Vector3[numOutputVerts];
            for (int i = 0; i < numOutputVerts; i++)
            {
                vertices[i] = outputVerts[i];
            }

            meshRef.Mesh.vertices = vertices;
            meshRef.Mesh.RecalculateNormals();
            meshRef.Mesh.RecalculateBounds();

            outputVerts.Dispose();
        }

        entities.Dispose();
    }
}

/// <summary>
/// Cloth Mesh更新系统 - 在Presentation阶段运行
/// </summary>
[UpdateInGroup(typeof(PresentationSystemGroup))]
public partial class ClothMeshUpdateSystem : SystemBase
{
    private EntityQuery _clothQuery;

    protected override void OnCreate()
    {
        _clothQuery = GetEntityQuery(
            ComponentType.ReadOnly<ClothTag>(),
            ComponentType.ReadOnly<ClothSolverConfig>(),
            ComponentType.ReadOnly<ParticlePosition>()
        );
        RequireForUpdate(_clothQuery);
    }

    protected override void OnUpdate()
    {
        var entities = _clothQuery.ToEntityArray(Allocator.Temp);

        for (int e = 0; e < entities.Length; e++)
        {
            var entity = entities[e];

            if (!EntityManager.HasComponent<ManagedMeshReference>(entity)) continue;
            var meshRef = EntityManager.GetComponentObject<ManagedMeshReference>(entity);
            if (meshRef == null || meshRef.Mesh == null) continue;

            var cfg = EntityManager.GetComponentData<ClothSolverConfig>(entity);
            var positions = EntityManager.GetBuffer<ParticlePosition>(entity);

            // 数据有效性检查
            if (positions.Length == 0 || cfg.NumParticles <= 0) continue;

            var posArr = positions.Reinterpret<float3>().AsNativeArray();
            int numParticles = cfg.NumParticles;

            var vertices = new Vector3[numParticles];
            for (int i = 0; i < numParticles; i++)
            {
                vertices[i] = posArr[i];
            }

            meshRef.Mesh.vertices = vertices;
            meshRef.Mesh.RecalculateNormals();
            meshRef.Mesh.RecalculateTangents();
            meshRef.Mesh.RecalculateBounds();
        }

        entities.Dispose();
    }
}

/// <summary>
/// 托管Mesh引用组件（IComponentData不能存储托管类型，使用IComponentData的class版本）
/// </summary>
public class ManagedMeshReference : IComponentData
{
    public Mesh Mesh;
    public MeshFilter MeshFilter;
    public MeshRenderer MeshRenderer;
}
