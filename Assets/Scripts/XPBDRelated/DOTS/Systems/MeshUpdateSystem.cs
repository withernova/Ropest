using Unity.Burst;
using Unity.Collections;
using Unity.Entities;
using Unity.Jobs;
using Unity.Mathematics;
using UnityEngine;
using UnityEngine.Rendering;

[BurstCompile]
public struct ApplyWorldToLocalPositionsJob : IJobParallelFor
{
    public NativeArray<float3> Positions;
    public float4x4 WorldToLocal;

    public void Execute(int i)
    {
        Positions[i] = math.transform(WorldToLocal, Positions[i]);
    }
}

[BurstCompile]
public struct ApplyWorldToLocalNormalsJob : IJobParallelFor
{
    public NativeArray<float3> Normals;
    public float3x3 WorldToLocal3x3;

    public void Execute(int i)
    {
        float3 n = math.mul(WorldToLocal3x3, Normals[i]);
        float len = math.length(n);
        Normals[i] = len > 1e-8f ? n / len : Normals[i];
    }
}

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
        _ropeQuery.CompleteDependency();

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

            // 粒子位置是 world 空间，写回 Mesh 前必须转换到 MeshFilter.transform 的 local 空间，
            // 否则当 GameObject 非 identity 时视觉位置与仿真/碰撞位置错位。
            if (meshRef.MeshFilter != null)
            {
                var worldToLocal = (float4x4)meshRef.MeshFilter.transform.worldToLocalMatrix;
                new ApplyWorldToLocalPositionsJob
                {
                    Positions = _outputVerts,
                    WorldToLocal = worldToLocal
                }.Schedule(numOutputVerts, 64).Complete();
            }

            // 使用 NativeArray 版 SetVertices 直接把 float3 作为 Vector3 传入，零托管分配
            var vertsAsVec3 = _outputVerts.Reinterpret<Vector3>();
            meshRef.Mesh.SetVertices(vertsAsVec3);
            meshRef.Mesh.RecalculateNormals();
            meshRef.Mesh.RecalculateBounds();
        }

        entities.Dispose();
    }
}

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
        _clothQuery.CompleteDependency();

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

                // 粒子位置/法线都在 world 空间，需转到 MeshFilter local 空间写回 Mesh。
                if (meshRef.MeshFilter != null)
                {
                    var worldToLocal = (float4x4)meshRef.MeshFilter.transform.worldToLocalMatrix;
                    var w2l3x3 = new float3x3(worldToLocal);
                    var posHandle = new ApplyWorldToLocalPositionsJob
                    {
                        Positions = _renderPositions,
                        WorldToLocal = worldToLocal
                    }.Schedule(numRenderVerts, 64);
                    var nrmHandle = new ApplyWorldToLocalNormalsJob
                    {
                        Normals = _renderNormals,
                        WorldToLocal3x3 = w2l3x3
                    }.Schedule(numRenderVerts, 64);
                    JobHandle.CombineDependencies(posHandle, nrmHandle).Complete();
                }

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
                EnsureCapacity(ref _renderPositions, numParticles);
                NativeArray<float3>.Copy(posArr, 0, _renderPositions, 0, numParticles);

                if (meshRef.MeshFilter != null)
                {
                    var worldToLocal = (float4x4)meshRef.MeshFilter.transform.worldToLocalMatrix;
                    new ApplyWorldToLocalPositionsJob
                    {
                        Positions = _renderPositions,
                        WorldToLocal = worldToLocal
                    }.Schedule(numParticles, 64).Complete();
                }

                var posVec = _renderPositions.Reinterpret<Vector3>();
                // SetVertices 第 3 个参数 length 用 numParticles，确保 _renderPositions 容量更大时不越界
                meshRef.Mesh.SetVertices(posVec, 0, numParticles);
                meshRef.Mesh.RecalculateNormals();
                meshRef.Mesh.RecalculateBounds();
            }
        }

        entities.Dispose();
    }
}

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
        _softBodyQuery.CompleteDependency();

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

                // world → local 变换
                if (meshRef.MeshFilter != null)
                {
                    var worldToLocal = (float4x4)meshRef.MeshFilter.transform.worldToLocalMatrix;
                    var w2l3x3 = new float3x3(worldToLocal);
                    var posHandle = new ApplyWorldToLocalPositionsJob
                    {
                        Positions = _renderPositions,
                        WorldToLocal = worldToLocal
                    }.Schedule(numRenderVerts, 64);
                    var nrmHandle = new ApplyWorldToLocalNormalsJob
                    {
                        Normals = _renderNormals,
                        WorldToLocal3x3 = w2l3x3
                    }.Schedule(numRenderVerts, 64);
                    JobHandle.CombineDependencies(posHandle, nrmHandle).Complete();
                }

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
                EnsureCapacity(ref _renderPositions, numParticles);
                NativeArray<float3>.Copy(posArr, 0, _renderPositions, 0, numParticles);

                if (meshRef.MeshFilter != null)
                {
                    var worldToLocal = (float4x4)meshRef.MeshFilter.transform.worldToLocalMatrix;
                    new ApplyWorldToLocalPositionsJob
                    {
                        Positions = _renderPositions,
                        WorldToLocal = worldToLocal
                    }.Schedule(numParticles, 64).Complete();
                }

                var posVec = _renderPositions.Reinterpret<Vector3>();
                meshRef.Mesh.SetVertices(posVec, 0, numParticles);
                meshRef.Mesh.RecalculateNormals();
                meshRef.Mesh.RecalculateBounds();
            }
        }

        entities.Dispose();
    }
}

public class ManagedMeshReference : IComponentData
{
    public Mesh Mesh;
    public MeshFilter MeshFilter;
    public MeshRenderer MeshRenderer;
}
