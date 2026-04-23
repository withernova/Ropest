using Unity.Collections;
using Unity.Entities;
using Unity.Jobs;
using Unity.Mathematics;
using UnityEngine;
using UnityEngine.Rendering;

/// <summary>
/// Rope Mesh更新系统 - 在Presentation阶段运行
/// 将ECS中的粒子数据回写到Unity Mesh
/// </summary>
[UpdateInGroup(typeof(PresentationSystemGroup))]
public partial class RopeMeshUpdateSystem : SystemBase
{
    private EntityQuery _ropeQuery;

    // 复用缓冲，避免每帧 new NativeArray / Vector3[]
    private NativeArray<float3> _outputVerts;

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

    protected override void OnDestroy()
    {
        if (_outputVerts.IsCreated) _outputVerts.Dispose();
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

            // 复用 NativeArray（容量不足时才重建）
            if (!_outputVerts.IsCreated || _outputVerts.Length != numOutputVerts)
            {
                if (_outputVerts.IsCreated) _outputVerts.Dispose();
                _outputVerts = new NativeArray<float3>(numOutputVerts, Allocator.Persistent);
            }

            var renderJob = new RopeRenderingJob
            {
                PointPos = posArr,
                GhostPos = ghostArr,
                OutputVertices = _outputVerts,
                Subdivision = subdivision,
                Radius = radius,
                NumPoints = numPoints
            };
            renderJob.Schedule(default(JobHandle)).Complete();

            // 使用 NativeArray 版 SetVertices 直接把 float3 作为 Vector3 传入，零托管分配
            var vertsAsVec3 = _outputVerts.Reinterpret<Vector3>();
            meshRef.Mesh.SetVertices(vertsAsVec3);
            meshRef.Mesh.RecalculateNormals();
            meshRef.Mesh.RecalculateBounds();
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

    // 复用缓冲，避免每帧 new
    private NativeArray<float3> _simNormals;
    private NativeArray<float3> _renderPositions;
    private NativeArray<float3> _renderNormals;

    protected override void OnCreate()
    {
        _clothQuery = GetEntityQuery(
            ComponentType.ReadOnly<ClothTag>(),
            ComponentType.ReadOnly<ClothSolverConfig>(),
            ComponentType.ReadOnly<ParticlePosition>()
        );
        RequireForUpdate(_clothQuery);
    }

    protected override void OnDestroy()
    {
        if (_simNormals.IsCreated) _simNormals.Dispose();
        if (_renderPositions.IsCreated) _renderPositions.Dispose();
        if (_renderNormals.IsCreated) _renderNormals.Dispose();
    }

    private static void EnsureCapacity(ref NativeArray<float3> arr, int size)
    {
        if (!arr.IsCreated || arr.Length != size)
        {
            if (arr.IsCreated) arr.Dispose();
            arr = new NativeArray<float3>(size, Allocator.Persistent);
        }
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

                int numParticles = cfg.NumParticles;
                EnsureCapacity(ref _simNormals, numParticles);
                EnsureCapacity(ref _renderPositions, numRenderVerts);
                EnsureCapacity(ref _renderNormals, numRenderVerts);

                // 第一步：计算模拟网格的平滑顶点法线
                var normalJob = new ComputeSimNormalsJob
                {
                    Positions = posArr,
                    Triangles = triArr,
                    Normals = _simNormals
                };
                var normalHandle = normalJob.Schedule();

                // 绑定数据：零拷贝复用 Buffer 内存（布局与 RenderVertexBinding 一致）
                var bindingArr = bindingBuf.Reinterpret<RenderBindingData>().AsNativeArray();

                // 第二步：同时插值位置和法线（依赖 normalJob 完成）
                var interpJob = new InterpolatePositionAndNormalJob
                {
                    SimPositions = posArr,
                    SimNormals = _simNormals,
                    Bindings = bindingArr,
                    RenderPositions = _renderPositions,
                    RenderNormals = _renderNormals
                };
                interpJob.Schedule(numRenderVerts, 64, normalHandle).Complete();

                // 回写：使用 NativeArray 版本的 SetVertices/SetNormals，避免托管数组分配和拷贝
                var renderPosVec = _renderPositions.Reinterpret<Vector3>();
                var renderNormVec = _renderNormals.Reinterpret<Vector3>();
                meshRef.Mesh.SetVertices(renderPosVec);
                meshRef.Mesh.SetNormals(renderNormVec);
                meshRef.Mesh.RecalculateBounds();
            }
            else
            {
                int numParticles = cfg.NumParticles;
                // 直接用粒子位置作为顶点，零拷贝 Reinterpret 传入
                var posVec = posArr.Reinterpret<Vector3>();
                meshRef.Mesh.SetVertices(posVec);
                meshRef.Mesh.RecalculateNormals();
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

    // 复用缓冲，避免每帧 new
    private NativeArray<float3> _simNormals;
    private NativeArray<float3> _renderPositions;
    private NativeArray<float3> _renderNormals;

    protected override void OnCreate()
    {
        _softBodyQuery = GetEntityQuery(
            ComponentType.ReadOnly<SoftBodyTag>(),
            ComponentType.ReadOnly<SoftBodySolverConfig>(),
            ComponentType.ReadOnly<ParticlePosition>()
        );
        RequireForUpdate(_softBodyQuery);
    }

    protected override void OnDestroy()
    {
        if (_simNormals.IsCreated) _simNormals.Dispose();
        if (_renderPositions.IsCreated) _renderPositions.Dispose();
        if (_renderNormals.IsCreated) _renderNormals.Dispose();
    }

    private static void EnsureCapacity(ref NativeArray<float3> arr, int size)
    {
        if (!arr.IsCreated || arr.Length != size)
        {
            if (arr.IsCreated) arr.Dispose();
            arr = new NativeArray<float3>(size, Allocator.Persistent);
        }
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

                int numParticles = cfg.NumParticles;
                EnsureCapacity(ref _simNormals, numParticles);
                EnsureCapacity(ref _renderPositions, numRenderVerts);
                EnsureCapacity(ref _renderNormals, numRenderVerts);

                // 第一步：计算模拟网格的平滑顶点法线
                var normalJob = new ComputeSimNormalsJob
                {
                    Positions = posArr,
                    Triangles = triArr,
                    Normals = _simNormals
                };
                var normalHandle = normalJob.Schedule();

                // 绑定数据：零拷贝复用 Buffer 内存
                var bindingArr = bindingBuf.Reinterpret<RenderBindingData>().AsNativeArray();

                // 第二步：同时插值位置和法线（依赖 normalJob 完成）
                var interpJob = new InterpolatePositionAndNormalJob
                {
                    SimPositions = posArr,
                    SimNormals = _simNormals,
                    Bindings = bindingArr,
                    RenderPositions = _renderPositions,
                    RenderNormals = _renderNormals
                };
                interpJob.Schedule(numRenderVerts, 64, normalHandle).Complete();

                // 回写：使用 NativeArray 版本的 SetVertices/SetNormals，避免托管数组分配和拷贝
                var renderPosVec = _renderPositions.Reinterpret<Vector3>();
                var renderNormVec = _renderNormals.Reinterpret<Vector3>();
                meshRef.Mesh.SetVertices(renderPosVec);
                meshRef.Mesh.SetNormals(renderNormVec);
                meshRef.Mesh.RecalculateBounds();
            }
            else
            {
                int numParticles = cfg.NumParticles;
                // 直接用粒子位置作为顶点，零拷贝 Reinterpret 传入
                var posVec = posArr.Reinterpret<Vector3>();
                meshRef.Mesh.SetVertices(posVec);
                meshRef.Mesh.RecalculateNormals();
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
