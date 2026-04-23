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

            if (!EntityManager.HasComponent<ManagedMeshReference>(entity)) continue;
            var meshRef = EntityManager.GetComponentObject<ManagedMeshReference>(entity);
            if (meshRef == null || meshRef.Mesh == null) continue;

            var cfg = EntityManager.GetComponentData<RopeSolverConfig>(entity);
            var positions = EntityManager.GetBuffer<ParticlePosition>(entity);
            var ghostPositions = EntityManager.GetBuffer<GhostPosition>(entity);

            if (positions.Length == 0 || ghostPositions.Length == 0 || cfg.NumPoints <= 0) continue;

            var posArr = positions.Reinterpret<float3>().AsNativeArray();
            var ghostArr = ghostPositions.Reinterpret<float3>().AsNativeArray();

            int numPoints = cfg.NumPoints;
            int subdivision = cfg.Subdivision;
            float radius = cfg.Radius;

            int interpolations = 10;
            int numOutputVerts = 1 + (numPoints - 2) * (interpolations - 1) * subdivision + 1;

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

            var vertices = new Vector3[numOutputVerts];
            for (int i = 0; i < numOutputVerts; i++)
                vertices[i] = outputVerts[i];

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
/// 支持两种模式：
/// 1. 直接模式：模拟粒子位置直接写回Mesh顶点
/// 2. 细分平滑模式：位置线性插值 + 法线平滑插值，实现视觉平滑
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

            if (positions.Length == 0 || cfg.NumParticles <= 0) continue;

            var posArr = positions.Reinterpret<float3>().AsNativeArray();

            bool useSubdivision = EntityManager.HasComponent<RenderMeshConfig>(entity) &&
                                  EntityManager.GetComponentData<RenderMeshConfig>(entity).UseSubdivision &&
                                  EntityManager.HasComponent<RenderVertexBinding>(entity);

            if (useSubdivision)
            {
                var renderCfg = EntityManager.GetComponentData<RenderMeshConfig>(entity);
                var bindingBuf = EntityManager.GetBuffer<RenderVertexBinding>(entity);
                int numRenderVerts = renderCfg.NumRenderVertices;

                if (bindingBuf.Length == 0 || numRenderVerts <= 0) continue;

                // 获取模拟三角形数据（零拷贝复用 Buffer 内存）
                var triBuf = EntityManager.GetBuffer<TriangleIndex>(entity);
                if (triBuf.Length == 0) continue;
                var triArr = triBuf.Reinterpret<int3>().AsNativeArray();

                // 第一步：计算模拟网格的平滑顶点法线
                int numParticles = cfg.NumParticles;
                var simNormals = new NativeArray<float3>(numParticles, Allocator.TempJob);
                var normalJob = new ComputeSimNormalsJob
                {
                    Positions = posArr,
                    Triangles = triArr,
                    Normals = simNormals
                };
                var normalHandle = normalJob.Schedule();

                // 绑定数据：零拷贝复用 Buffer 内存（布局与 RenderVertexBinding 一致）
                var bindingArr = bindingBuf.Reinterpret<RenderBindingData>().AsNativeArray();

                var renderPositions = new NativeArray<float3>(numRenderVerts, Allocator.TempJob);
                var renderNormals = new NativeArray<float3>(numRenderVerts, Allocator.TempJob);

                // 第二步：同时插值位置和法线（依赖 normalJob 完成）
                var interpJob = new InterpolatePositionAndNormalJob
                {
                    SimPositions = posArr,
                    SimNormals = simNormals,
                    Bindings = bindingArr,
                    RenderPositions = renderPositions,
                    RenderNormals = renderNormals
                };
                interpJob.Schedule(numRenderVerts, 64, normalHandle).Complete();

                // 回写：使用 NativeArray 版本的 SetVertices/SetNormals，避免托管数组分配和拷贝
                var renderPosVec = renderPositions.Reinterpret<Vector3>();
                var renderNormVec = renderNormals.Reinterpret<Vector3>();
                meshRef.Mesh.SetVertices(renderPosVec);
                meshRef.Mesh.SetNormals(renderNormVec);
                meshRef.Mesh.RecalculateBounds();

                renderPositions.Dispose();
                renderNormals.Dispose();
                simNormals.Dispose();
            }
            else
            {
                int numParticles = cfg.NumParticles;
                var vertices = new Vector3[numParticles];
                for (int i = 0; i < numParticles; i++)
                    vertices[i] = posArr[i];

                meshRef.Mesh.vertices = vertices;
                meshRef.Mesh.RecalculateNormals();
                meshRef.Mesh.RecalculateTangents();
                meshRef.Mesh.RecalculateBounds();
            }
        }

        entities.Dispose();
    }
}

/// <summary>
/// SoftBody Mesh更新系统 - 在Presentation阶段运行
/// 支持两种模式：
/// 1. 直接模式：模拟粒子位置直接写回Mesh顶点
/// 2. 细分平滑模式：位置线性插值 + 法线平滑插值，实现视觉平滑
/// </summary>
[UpdateInGroup(typeof(PresentationSystemGroup))]
public partial class SoftBodyMeshUpdateSystem : SystemBase
{
    private EntityQuery _softBodyQuery;

    protected override void OnCreate()
    {
        _softBodyQuery = GetEntityQuery(
            ComponentType.ReadOnly<SoftBodyTag>(),
            ComponentType.ReadOnly<SoftBodySolverConfig>(),
            ComponentType.ReadOnly<ParticlePosition>()
        );
        RequireForUpdate(_softBodyQuery);
    }

    protected override void OnUpdate()
    {
        var entities = _softBodyQuery.ToEntityArray(Allocator.Temp);

        for (int e = 0; e < entities.Length; e++)
        {
            var entity = entities[e];

            if (!EntityManager.HasComponent<ManagedMeshReference>(entity)) continue;
            var meshRef = EntityManager.GetComponentObject<ManagedMeshReference>(entity);
            if (meshRef == null || meshRef.Mesh == null) continue;

            var cfg = EntityManager.GetComponentData<SoftBodySolverConfig>(entity);
            var positions = EntityManager.GetBuffer<ParticlePosition>(entity);

            if (positions.Length == 0 || cfg.NumParticles <= 0) continue;

            var posArr = positions.Reinterpret<float3>().AsNativeArray();

            bool useSubdivision = EntityManager.HasComponent<RenderMeshConfig>(entity) &&
                                  EntityManager.GetComponentData<RenderMeshConfig>(entity).UseSubdivision &&
                                  EntityManager.HasComponent<RenderVertexBinding>(entity);

            if (useSubdivision)
            {
                var renderCfg = EntityManager.GetComponentData<RenderMeshConfig>(entity);
                var bindingBuf = EntityManager.GetBuffer<RenderVertexBinding>(entity);
                int numRenderVerts = renderCfg.NumRenderVertices;

                if (bindingBuf.Length == 0 || numRenderVerts <= 0) continue;

                // 获取表面三角形数据（零拷贝复用 Buffer 内存）
                var surfTriBuf = EntityManager.GetBuffer<SurfaceTriangleIndex>(entity);
                if (surfTriBuf.Length == 0) continue;
                var triArr = surfTriBuf.Reinterpret<int3>().AsNativeArray();

                // 第一步：计算模拟网格的平滑顶点法线
                int numParticles = cfg.NumParticles;
                var simNormals = new NativeArray<float3>(numParticles, Allocator.TempJob);
                var normalJob = new ComputeSimNormalsJob
                {
                    Positions = posArr,
                    Triangles = triArr,
                    Normals = simNormals
                };
                var normalHandle = normalJob.Schedule();

                // 绑定数据：零拷贝复用 Buffer 内存
                var bindingArr = bindingBuf.Reinterpret<RenderBindingData>().AsNativeArray();

                var renderPositions = new NativeArray<float3>(numRenderVerts, Allocator.TempJob);
                var renderNormals = new NativeArray<float3>(numRenderVerts, Allocator.TempJob);

                // 第二步：同时插值位置和法线（依赖 normalJob 完成）
                var interpJob = new InterpolatePositionAndNormalJob
                {
                    SimPositions = posArr,
                    SimNormals = simNormals,
                    Bindings = bindingArr,
                    RenderPositions = renderPositions,
                    RenderNormals = renderNormals
                };
                interpJob.Schedule(numRenderVerts, 64, normalHandle).Complete();

                // 回写：使用 NativeArray 版本的 SetVertices/SetNormals，避免托管数组分配和拷贝
                var renderPosVec = renderPositions.Reinterpret<Vector3>();
                var renderNormVec = renderNormals.Reinterpret<Vector3>();
                meshRef.Mesh.SetVertices(renderPosVec);
                meshRef.Mesh.SetNormals(renderNormVec);
                meshRef.Mesh.RecalculateBounds();

                renderPositions.Dispose();
                renderNormals.Dispose();
                simNormals.Dispose();
            }
            else
            {
                int numParticles = cfg.NumParticles;
                var vertices = new Vector3[numParticles];
                for (int i = 0; i < numParticles; i++)
                    vertices[i] = posArr[i];

                meshRef.Mesh.vertices = vertices;
                meshRef.Mesh.RecalculateNormals();
                meshRef.Mesh.RecalculateTangents();
                meshRef.Mesh.RecalculateBounds();
            }
        }

        entities.Dispose();
    }
}

/// <summary>
/// 托管Mesh引用组件
/// </summary>
public class ManagedMeshReference : IComponentData
{
    public Mesh Mesh;
    public MeshFilter MeshFilter;
    public MeshRenderer MeshRenderer;
}
