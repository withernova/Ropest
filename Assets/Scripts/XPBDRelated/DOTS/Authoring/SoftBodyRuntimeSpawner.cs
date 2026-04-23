using System.Collections.Generic;
using Unity.Collections;
using Unity.Entities;
using Unity.Mathematics;
using UnityEngine;

/// <summary>
/// 运行时SoftBody初始化器 - 不依赖SubScene，可以在运行时动态创建SoftBody Entity
/// 支持两种模式：
/// 1. 自动生成立方体软体（通过分辨率参数控制）
/// 2. 使用外部Mesh（需要提供四面体数据或自动进行简单四面体化）
/// 
/// 使用方式：
/// 1. 创建空GameObject，挂载此脚本
/// 2. 设置参数后运行，自动创建软体Entity
/// </summary>
public class SoftBodyRuntimeSpawner : MonoBehaviour
{
    [Header("环境参数")]
    public Vector3 gravity = new(0, -9.8f, 0);

    public enum SoftBodyShape
    {
        Cuboid,
        Sphere
    }

    [Header("软体形状")]
    [Tooltip("软体形状：立方体 / 球体（椭球）")]
    public SoftBodyShape shape = SoftBodyShape.Sphere;
    [Tooltip("X方向尺寸（球体时为 X 轴半径 * 2）")]
    [Range(0.1f, 20f)] public float sizeX = 1f;
    [Tooltip("Y方向尺寸（球体时为 Y 轴半径 * 2）")]
    [Range(0.1f, 20f)] public float sizeY = 1f;
    [Tooltip("Z方向尺寸（球体时为 Z 轴半径 * 2）")]
    [Range(0.1f, 20f)] public float sizeZ = 1f;
    [Tooltip("每个轴方向的分段数")]
    [Range(1, 20)] public int resolution = 4;
    [Tooltip("球体模式下，是否把表面顶点投影到椭球面（更圆滑）")]
    public bool projectToSphereSurface = true;
    public Vector3 meshOrigin;
    public Material softBodyMaterial;

    [Header("求解器参数")]
    [Range(1, 30)] public int numSubSteps = 10;
    [Tooltip("距离约束柔度（0=完全刚性，越大越软）")]
    [Range(0f, 1f)] public float distanceStiffness = 0f;
    [Tooltip("体积约束柔度（0=不可压缩，越大越容易压缩）")]
    [Range(0f, 1f)] public float volumeStiffness = 0f;
    [Range(0f, 1f)] public float damping = 0.02f;
    [Range(0.01f, 0.5f)] public float collisionRadius = 0.05f;
    [Range(0f, 1f)] public float friction = 0.3f;

    [Header("渲染网格细分")]
    [Tooltip("渲染网格细分迭代次数（0=不细分直接用模拟网格，1=4倍面数，2=16倍面数）")]
    [Range(0, 3)] public int renderSubdivisionIterations = 1;

    [Header("固定点")]
    [Tooltip("固定的顶点索引列表（invMass设为0）")]
    public List<int> fixedVertices = new List<int>();

    [Header("固定顶部")]
    [Tooltip("是否固定Y坐标最大的一层顶点")]
    public bool fixTopLayer = false;

    private Entity softBodyEntity;
    private EntityManager entityManager;
    private Mesh surfaceMesh;
    private bool entityCreated = false;

    void Start()
    {
        entityManager = World.DefaultGameObjectInjectionWorld.EntityManager;
        SpawnSoftBody();
    }

    void SpawnSoftBody()
    {
        // 根据形状生成网格数据（顶点 + 四面体 + 表面三角形 + 边）
        Vector3[] vertices;
        int[] tetrahedra;       // 每4个为一组
        int[] surfaceTriangles; // 每3个为一组
        int[] edges;            // 每2个为一组（去重）

        if (shape == SoftBodyShape.Sphere)
        {
            GenerateSphereMeshData(out vertices, out tetrahedra, out surfaceTriangles, out edges);
        }
        else
        {
            GenerateCuboidMeshData(out vertices, out tetrahedra, out surfaceTriangles, out edges);
        }

        int numParticles = vertices.Length;
        int numTets = tetrahedra.Length / 4;
        int numEdges = edges.Length / 2;

        // 创建表面Mesh用于渲染（支持细分高面数渲染）
        bool useSubdivision = renderSubdivisionIterations > 0;
        Vector3[] renderVertices;
        int[] renderTriangles;
        Vector2[] renderUVs;
        SubdivisionMeshHelper.BindingData[] bindingData = null;

        if (useSubdivision)
        {
            // 提取表面顶点和UV用于细分
            var surfVerts = new Vector3[numParticles];
            System.Array.Copy(vertices, surfVerts, numParticles);
            var surfUVs = GenerateSphereUVs(surfVerts);

            // 细分的同时精确追踪每个顶点的重心坐标绑定（无搜索误差）
            SubdivisionMeshHelper.SubdivideWithBindings(
                surfVerts, surfaceTriangles, surfUVs, renderSubdivisionIterations,
                out renderVertices, out renderTriangles, out renderUVs, out bindingData);

            Debug.Log($"[SoftBody] 细分渲染网格: 模拟顶点={numParticles}, 渲染顶点={renderVertices.Length}, " +
                      $"模拟三角形={surfaceTriangles.Length / 3}, 渲染三角形={renderTriangles.Length / 3}");
        }
        else
        {
            renderVertices = vertices;
            renderTriangles = surfaceTriangles;
            renderUVs = GenerateSphereUVs(vertices);
        }

        surfaceMesh = CreateSurfaceMeshFromData(renderVertices, renderTriangles, renderUVs);

        var meshFilter = gameObject.AddComponent<MeshFilter>();
        var meshRenderer = gameObject.AddComponent<MeshRenderer>();
        meshFilter.mesh = surfaceMesh;
        meshRenderer.material = softBodyMaterial;

        // 计算逆质量（基于相邻四面体体积）
        var invMasses = new float[numParticles];
        for (int i = 0; i < numParticles; i++) invMasses[i] = 0f;

        // 累加每个顶点关联的四面体体积
        for (int t = 0; t < numTets; t++)
        {
            int i0 = tetrahedra[t * 4 + 0];
            int i1 = tetrahedra[t * 4 + 1];
            int i2 = tetrahedra[t * 4 + 2];
            int i3 = tetrahedra[t * 4 + 3];

            float vol = ComputeTetVolume(vertices[i0], vertices[i1], vertices[i2], vertices[i3]);
            float absVol = Mathf.Abs(vol);

            invMasses[i0] += absVol / 4f;
            invMasses[i1] += absVol / 4f;
            invMasses[i2] += absVol / 4f;
            invMasses[i3] += absVol / 4f;
        }

        for (int i = 0; i < numParticles; i++)
        {
            invMasses[i] = invMasses[i] > 1e-8f ? 1f / invMasses[i] : 1f;
        }

        // 固定指定顶点
        foreach (int idx in fixedVertices)
        {
            if (idx >= 0 && idx < numParticles)
                invMasses[idx] = 0f;
        }

        // 固定顶部一层
        if (fixTopLayer)
        {
            float maxY = float.MinValue;
            for (int i = 0; i < numParticles; i++)
                maxY = Mathf.Max(maxY, vertices[i].y);

            float threshold = maxY - (sizeY / resolution) * 0.5f;
            for (int i = 0; i < numParticles; i++)
            {
                if (vertices[i].y >= threshold)
                    invMasses[i] = 0f;
            }
        }

        // 计算四面体静止体积
        var restVolumes = new float[numTets];
        for (int t = 0; t < numTets; t++)
        {
            int i0 = tetrahedra[t * 4 + 0];
            int i1 = tetrahedra[t * 4 + 1];
            int i2 = tetrahedra[t * 4 + 2];
            int i3 = tetrahedra[t * 4 + 3];
            restVolumes[t] = ComputeTetVolume(vertices[i0], vertices[i1], vertices[i2], vertices[i3]);
        }

        // 计算边的静止长度
        var restLengths = new float[numEdges];
        for (int e2 = 0; e2 < numEdges; e2++)
        {
            int a = edges[e2 * 2];
            int b = edges[e2 * 2 + 1];
            restLengths[e2] = Vector3.Distance(vertices[a], vertices[b]);
        }

        // === 使用Archetype一次性创建Entity ===
        var componentTypes = new List<ComponentType>
        {
            typeof(SoftBodyTag),
            typeof(MeshUpdateTag),
            typeof(SoftBodySolverConfig),
            typeof(ParticlePosition),
            typeof(ParticlePrevPosition),
            typeof(ParticleVelocity),
            typeof(ParticleInvMass),
            typeof(SoftBodyEdge),
            typeof(SoftBodyDistanceLambda),
            typeof(Tetrahedron),
            typeof(TetrahedronRestVolume),
            typeof(TetrahedronVolumeLambda),
            typeof(SurfaceTriangleIndex),
            typeof(RenderMeshConfig)
        };
        if (useSubdivision)
        {
            componentTypes.Add(typeof(RenderVertexBinding));
        }
        var archetype = entityManager.CreateArchetype(componentTypes.ToArray());

        softBodyEntity = entityManager.CreateEntity(archetype);

        // 设置配置
        entityManager.SetComponentData(softBodyEntity, new SoftBodySolverConfig
        {
            NumParticles = numParticles,
            NumSubSteps = numSubSteps,
            Gravity = gravity,
            DistanceStiffness = distanceStiffness,
            VolumeStiffness = volumeStiffness,
            Damping = damping,
            CollisionRadius = collisionRadius,
            Friction = friction
        });

        // === 填充Buffer数据 ===
        var posBuf = entityManager.GetBuffer<ParticlePosition>(softBodyEntity);
        var prevBuf = entityManager.GetBuffer<ParticlePrevPosition>(softBodyEntity);
        var velBuf = entityManager.GetBuffer<ParticleVelocity>(softBodyEntity);
        var massBuf = entityManager.GetBuffer<ParticleInvMass>(softBodyEntity);

        for (int i = 0; i < numParticles; i++)
        {
            float3 pos = vertices[i];
            posBuf.Add(new ParticlePosition { Value = pos });
            prevBuf.Add(new ParticlePrevPosition { Value = pos });
            velBuf.Add(new ParticleVelocity { Value = float3.zero });
            massBuf.Add(new ParticleInvMass { Value = invMasses[i] });
        }

        // 边数据
        var edgeBuf = entityManager.GetBuffer<SoftBodyEdge>(softBodyEntity);
        var distLambdaBuf = entityManager.GetBuffer<SoftBodyDistanceLambda>(softBodyEntity);

        for (int e2 = 0; e2 < numEdges; e2++)
        {
            edgeBuf.Add(new SoftBodyEdge
            {
                IndexA = edges[e2 * 2],
                IndexB = edges[e2 * 2 + 1],
                RestLength = restLengths[e2]
            });
            distLambdaBuf.Add(new SoftBodyDistanceLambda { Value = 0f });
        }

        // 四面体数据
        var tetBuf = entityManager.GetBuffer<Tetrahedron>(softBodyEntity);
        var tetVolBuf = entityManager.GetBuffer<TetrahedronRestVolume>(softBodyEntity);
        var volLambdaBuf = entityManager.GetBuffer<TetrahedronVolumeLambda>(softBodyEntity);

        for (int t = 0; t < numTets; t++)
        {
            tetBuf.Add(new Tetrahedron
            {
                I0 = tetrahedra[t * 4 + 0],
                I1 = tetrahedra[t * 4 + 1],
                I2 = tetrahedra[t * 4 + 2],
                I3 = tetrahedra[t * 4 + 3]
            });
            tetVolBuf.Add(new TetrahedronRestVolume { Value = restVolumes[t] });
            volLambdaBuf.Add(new TetrahedronVolumeLambda { Value = 0f });
        }

        // 表面三角形索引（用于渲染插值时查找模拟三角形）
        var surfTriBuf = entityManager.GetBuffer<SurfaceTriangleIndex>(softBodyEntity);
        for (int i = 0; i < surfaceTriangles.Length; i += 3)
        {
            surfTriBuf.Add(new SurfaceTriangleIndex
            {
                I0 = surfaceTriangles[i],
                I1 = surfaceTriangles[i + 1],
                I2 = surfaceTriangles[i + 2]
            });
        }

        // 设置渲染网格配置
        entityManager.SetComponentData(softBodyEntity, new RenderMeshConfig
        {
            NumRenderVertices = renderVertices.Length,
            UseSubdivision = useSubdivision
        });

        // 填充渲染顶点绑定数据
        if (useSubdivision && bindingData != null)
        {
            var bindBuf = entityManager.GetBuffer<RenderVertexBinding>(softBodyEntity);
            for (int i = 0; i < bindingData.Length; i++)
            {
                var bd = bindingData[i];
                bindBuf.Add(new RenderVertexBinding
                {
                    SimI0 = bd.SimI0,
                    SimI1 = bd.SimI1,
                    SimI2 = bd.SimI2,
                    U = bd.U,
                    V = bd.V,
                    W = bd.W
                });
            }
        }

        // 托管Mesh引用
        entityManager.AddComponentObject(softBodyEntity, new ManagedMeshReference
        {
            Mesh = surfaceMesh,
            MeshFilter = GetComponent<MeshFilter>(),
            MeshRenderer = GetComponent<MeshRenderer>()
        });

        entityCreated = true;
    }

    void OnDestroy()
    {
        if (!entityCreated) return;
        if (World.DefaultGameObjectInjectionWorld == null || !World.DefaultGameObjectInjectionWorld.IsCreated) return;
        if (entityManager.Exists(softBodyEntity))
        {
            entityManager.DestroyEntity(softBodyEntity);
        }
    }

    /// <summary>
    /// 生成球体（椭球）网格数据：顶点、四面体、表面三角形、边
    /// 基于体素化思路：先构建 resolution^3 的立方网格，只保留其中心落在椭球内的小立方体，
    /// 每个保留的小立方体被切分为 5 个四面体；随后把表面顶点投影到椭球面以获得圆滑外观。
    /// 这样可以完整复用 Cuboid 模式的拓扑、边、体积、表面提取流程。
    /// </summary>
    void GenerateSphereMeshData(out Vector3[] vertices, out int[] tetrahedra,
        out int[] surfaceTriangles, out int[] edgesOut)
    {
        int nx = resolution;
        int ny = resolution;
        int nz = resolution;

        float dx = sizeX / nx;
        float dy = sizeY / ny;
        float dz = sizeZ / nz;

        // 椭球半径（以 meshOrigin 为一个角点，球心在 bbox 中心）
        float rx = sizeX * 0.5f;
        float ry = sizeY * 0.5f;
        float rz = sizeZ * 0.5f;
        Vector3 center = meshOrigin + new Vector3(rx, ry, rz);

        // 先生成全部候选网格顶点坐标（同立方体模式）
        int gridVertsCount = (nx + 1) * (ny + 1) * (nz + 1);
        var gridVerts = new Vector3[gridVertsCount];
        for (int ix = 0; ix <= nx; ix++)
        {
            for (int iy = 0; iy <= ny; iy++)
            {
                for (int iz = 0; iz <= nz; iz++)
                {
                    int idx = ix * (ny + 1) * (nz + 1) + iy * (nz + 1) + iz;
                    gridVerts[idx] = meshOrigin + new Vector3(ix * dx, iy * dy, iz * dz);
                }
            }
        }

        // 标记哪些小立方体的中心落在椭球内
        // 同时收集这些立方体用到的顶点
        bool IsInsideEllipsoid(Vector3 p)
        {
            Vector3 d = p - center;
            float fx = rx > 1e-6f ? d.x / rx : 0f;
            float fy = ry > 1e-6f ? d.y / ry : 0f;
            float fz = rz > 1e-6f ? d.z / rz : 0f;
            return (fx * fx + fy * fy + fz * fz) <= 1.0f + 1e-6f;
        }

        var usedVertexMap = new Dictionary<int, int>(); // oldIdx -> newIdx
        var newVerts = new List<Vector3>();
        var tetList = new List<int>();

        int MapVertex(int oldIdx)
        {
            if (!usedVertexMap.TryGetValue(oldIdx, out int newIdx))
            {
                newIdx = newVerts.Count;
                newVerts.Add(gridVerts[oldIdx]);
                usedVertexMap[oldIdx] = newIdx;
            }
            return newIdx;
        }

        // 如果分辨率过低导致没有任何小立方体的中心落入椭球（极小 resolution 时可能发生），
        // 则退化为保留所有 8 个角都在椭球内的立方体。若仍然没有，则至少保留 bbox 中心处的立方体。
        // 这里直接用“中心在内”即可。
        bool anyKept = false;
        for (int ix = 0; ix < nx; ix++)
        {
            for (int iy = 0; iy < ny; iy++)
            {
                for (int iz = 0; iz < nz; iz++)
                {
                    Vector3 cubeCenter = meshOrigin + new Vector3((ix + 0.5f) * dx, (iy + 0.5f) * dy, (iz + 0.5f) * dz);
                    if (!IsInsideEllipsoid(cubeCenter)) continue;
                    anyKept = true;

                    int v0 = ix * (ny + 1) * (nz + 1) + iy * (nz + 1) + iz;
                    int v1 = (ix + 1) * (ny + 1) * (nz + 1) + iy * (nz + 1) + iz;
                    int v2 = ix * (ny + 1) * (nz + 1) + iy * (nz + 1) + (iz + 1);
                    int v3 = (ix + 1) * (ny + 1) * (nz + 1) + iy * (nz + 1) + (iz + 1);
                    int v4 = ix * (ny + 1) * (nz + 1) + (iy + 1) * (nz + 1) + iz;
                    int v5 = (ix + 1) * (ny + 1) * (nz + 1) + (iy + 1) * (nz + 1) + iz;
                    int v6 = ix * (ny + 1) * (nz + 1) + (iy + 1) * (nz + 1) + (iz + 1);
                    int v7 = (ix + 1) * (ny + 1) * (nz + 1) + (iy + 1) * (nz + 1) + (iz + 1);

                    int n0 = MapVertex(v0);
                    int n1 = MapVertex(v1);
                    int n2 = MapVertex(v2);
                    int n3 = MapVertex(v3);
                    int n4 = MapVertex(v4);
                    int n5 = MapVertex(v5);
                    int n6 = MapVertex(v6);
                    int n7 = MapVertex(v7);

                    if ((ix + iy + iz) % 2 == 0)
                    {
                        tetList.AddRange(new[] { n0, n1, n3, n5 });
                        tetList.AddRange(new[] { n0, n3, n2, n6 });
                        tetList.AddRange(new[] { n0, n5, n4, n6 });
                        tetList.AddRange(new[] { n3, n5, n7, n6 });
                        tetList.AddRange(new[] { n0, n3, n5, n6 });
                    }
                    else
                    {
                        tetList.AddRange(new[] { n1, n0, n2, n4 });
                        tetList.AddRange(new[] { n1, n2, n3, n7 });
                        tetList.AddRange(new[] { n1, n4, n5, n7 });
                        tetList.AddRange(new[] { n2, n4, n6, n7 });
                        tetList.AddRange(new[] { n1, n2, n4, n7 });
                    }
                }
            }
        }

        // 如果没有任何立方体被保留（resolution 过小），回退到立方体模式
        if (!anyKept)
        {
            Debug.LogWarning("[SoftBody] 球体体素化结果为空（resolution 过小），回退到立方体生成。");
            GenerateCuboidMeshData(out vertices, out tetrahedra, out surfaceTriangles, out edgesOut);
            return;
        }

        vertices = newVerts.ToArray();
        tetrahedra = tetList.ToArray();
        int numTets = tetrahedra.Length / 4;

        // 提取所有边（去重）
        var edgeSet = new HashSet<(int, int)>();
        for (int t = 0; t < numTets; t++)
        {
            int i0 = tetrahedra[t * 4 + 0];
            int i1 = tetrahedra[t * 4 + 1];
            int i2 = tetrahedra[t * 4 + 2];
            int i3 = tetrahedra[t * 4 + 3];

            AddEdge(edgeSet, i0, i1);
            AddEdge(edgeSet, i0, i2);
            AddEdge(edgeSet, i0, i3);
            AddEdge(edgeSet, i1, i2);
            AddEdge(edgeSet, i1, i3);
            AddEdge(edgeSet, i2, i3);
        }

        var edgeList = new List<int>();
        foreach (var e in edgeSet)
        {
            edgeList.Add(e.Item1);
            edgeList.Add(e.Item2);
        }
        edgesOut = edgeList.ToArray();

        // 提取表面三角形（只被一个四面体引用的面）
        var faceCount = new Dictionary<(int, int, int), int>();
        var faceInfo = new Dictionary<(int, int, int), (int, int, int, int)>();

        for (int t = 0; t < numTets; t++)
        {
            int i0 = tetrahedra[t * 4 + 0];
            int i1 = tetrahedra[t * 4 + 1];
            int i2 = tetrahedra[t * 4 + 2];
            int i3 = tetrahedra[t * 4 + 3];

            AddFace(faceCount, faceInfo, i0, i1, i2, i3);
            AddFace(faceCount, faceInfo, i0, i1, i3, i2);
            AddFace(faceCount, faceInfo, i0, i2, i3, i1);
            AddFace(faceCount, faceInfo, i1, i2, i3, i0);
        }

        // 收集所有表面顶点（属于任一表面三角形的顶点）
        var surfaceVertexSet = new HashSet<int>();
        var surfTriList = new List<int>();
        foreach (var kvp in faceCount)
        {
            if (kvp.Value != 1) continue;

            var info = faceInfo[kvp.Key];
            int a = info.Item1;
            int b = info.Item2;
            int c = info.Item3;
            int opposite = info.Item4;

            Vector3 pa = vertices[a];
            Vector3 pb = vertices[b];
            Vector3 pc = vertices[c];
            Vector3 po = vertices[opposite];

            Vector3 faceNormal = Vector3.Cross(pb - pa, pc - pa);
            Vector3 toOpposite = po - pa;

            if (Vector3.Dot(faceNormal, toOpposite) > 0f)
            {
                surfTriList.Add(a);
                surfTriList.Add(c);
                surfTriList.Add(b);
            }
            else
            {
                surfTriList.Add(a);
                surfTriList.Add(b);
                surfTriList.Add(c);
            }

            surfaceVertexSet.Add(a);
            surfaceVertexSet.Add(b);
            surfaceVertexSet.Add(c);
        }

        surfaceTriangles = surfTriList.ToArray();

        // 把表面顶点投影到椭球面（保持内部顶点不动，避免四面体退化）
        if (projectToSphereSurface)
        {
            foreach (int vi in surfaceVertexSet)
            {
                Vector3 d = vertices[vi] - center;
                // 归一化到椭球表面：令 (d.x/rx)^2 + (d.y/ry)^2 + (d.z/rz)^2 = 1
                float fx = rx > 1e-6f ? d.x / rx : 0f;
                float fy = ry > 1e-6f ? d.y / ry : 0f;
                float fz = rz > 1e-6f ? d.z / rz : 0f;
                float denom = Mathf.Sqrt(fx * fx + fy * fy + fz * fz);
                if (denom < 1e-6f) continue;
                float scale = 1f / denom;
                Vector3 projected = center + new Vector3(d.x * scale, d.y * scale, d.z * scale);
                vertices[vi] = projected;
            }
        }
    }

    /// <summary>
    /// 生成立方体网格数据：顶点、四面体、表面三角形、边
    /// 每个小立方体被分割为5个四面体
    /// </summary>
    void GenerateCuboidMeshData(out Vector3[] vertices, out int[] tetrahedra,
        out int[] surfaceTriangles, out int[] edgesOut)
    {
        int nx = resolution;
        int ny = resolution;
        int nz = resolution;

        float dx = sizeX / nx;
        float dy = sizeY / ny;
        float dz = sizeZ / nz;

        // 生成顶点网格
        int numVerts = (nx + 1) * (ny + 1) * (nz + 1);
        vertices = new Vector3[numVerts];

        for (int ix = 0; ix <= nx; ix++)
        {
            for (int iy = 0; iy <= ny; iy++)
            {
                for (int iz = 0; iz <= nz; iz++)
                {
                    int idx = ix * (ny + 1) * (nz + 1) + iy * (nz + 1) + iz;
                    vertices[idx] = meshOrigin + new Vector3(ix * dx, iy * dy, iz * dz);
                }
            }
        }

        // 每个小立方体分割为5个四面体
        var tetList = new List<int>();

        for (int ix = 0; ix < nx; ix++)
        {
            for (int iy = 0; iy < ny; iy++)
            {
                for (int iz = 0; iz < nz; iz++)
                {
                    // 小立方体的8个顶点索引
                    //    6----7
                    //   /|   /|
                    //  4----5 |
                    //  | 2--|-3
                    //  |/   |/
                    //  0----1
                    int v0 = ix * (ny + 1) * (nz + 1) + iy * (nz + 1) + iz;
                    int v1 = (ix + 1) * (ny + 1) * (nz + 1) + iy * (nz + 1) + iz;
                    int v2 = ix * (ny + 1) * (nz + 1) + iy * (nz + 1) + (iz + 1);
                    int v3 = (ix + 1) * (ny + 1) * (nz + 1) + iy * (nz + 1) + (iz + 1);
                    int v4 = ix * (ny + 1) * (nz + 1) + (iy + 1) * (nz + 1) + iz;
                    int v5 = (ix + 1) * (ny + 1) * (nz + 1) + (iy + 1) * (nz + 1) + iz;
                    int v6 = ix * (ny + 1) * (nz + 1) + (iy + 1) * (nz + 1) + (iz + 1);
                    int v7 = (ix + 1) * (ny + 1) * (nz + 1) + (iy + 1) * (nz + 1) + (iz + 1);

                    // 使用交替方向的5-四面体分割（保证相邻立方体共享面）
                    // 根据 (ix+iy+iz) 的奇偶性选择不同的分割方式
                    if ((ix + iy + iz) % 2 == 0)
                    {
                        // 分割方式A：中心对角线 v0-v7
                        tetList.AddRange(new[] { v0, v1, v3, v5 });
                        tetList.AddRange(new[] { v0, v3, v2, v6 });
                        tetList.AddRange(new[] { v0, v5, v4, v6 });
                        tetList.AddRange(new[] { v3, v5, v7, v6 });
                        tetList.AddRange(new[] { v0, v3, v5, v6 });
                    }
                    else
                    {
                        // 分割方式B：中心对角线 v1-v6
                        tetList.AddRange(new[] { v1, v0, v2, v4 });
                        tetList.AddRange(new[] { v1, v2, v3, v7 });
                        tetList.AddRange(new[] { v1, v4, v5, v7 });
                        tetList.AddRange(new[] { v2, v4, v6, v7 });
                        tetList.AddRange(new[] { v1, v2, v4, v7 });
                    }
                }
            }
        }

        tetrahedra = tetList.ToArray();

        // 提取所有边（去重）
        var edgeSet = new HashSet<(int, int)>();
        int numTets = tetrahedra.Length / 4;
        for (int t = 0; t < numTets; t++)
        {
            int i0 = tetrahedra[t * 4 + 0];
            int i1 = tetrahedra[t * 4 + 1];
            int i2 = tetrahedra[t * 4 + 2];
            int i3 = tetrahedra[t * 4 + 3];

            AddEdge(edgeSet, i0, i1);
            AddEdge(edgeSet, i0, i2);
            AddEdge(edgeSet, i0, i3);
            AddEdge(edgeSet, i1, i2);
            AddEdge(edgeSet, i1, i3);
            AddEdge(edgeSet, i2, i3);
        }

        var edgeList = new List<int>();
        foreach (var edge in edgeSet)
        {
            edgeList.Add(edge.Item1);
            edgeList.Add(edge.Item2);
        }
        edgesOut = edgeList.ToArray();

        // 提取表面三角形
        // 表面三角形 = 只被一个四面体引用的三角形面
        // 使用排序后的三元组作为key来统计面的引用次数
        // 同时记录每个面所属四面体的对面顶点，用于后续修正法线方向
        var faceCount = new Dictionary<(int, int, int), int>();
        // 排序key -> (三角形三个顶点, 对面顶点索引)
        var faceInfo = new Dictionary<(int, int, int), (int, int, int, int)>();

        for (int t = 0; t < numTets; t++)
        {
            int i0 = tetrahedra[t * 4 + 0];
            int i1 = tetrahedra[t * 4 + 1];
            int i2 = tetrahedra[t * 4 + 2];
            int i3 = tetrahedra[t * 4 + 3];

            // 四面体的4个面，每个面记录对面的顶点用于法线修正
            AddFace(faceCount, faceInfo, i0, i1, i2, i3); // 面(0,1,2)，对面顶点3
            AddFace(faceCount, faceInfo, i0, i1, i3, i2); // 面(0,1,3)，对面顶点2
            AddFace(faceCount, faceInfo, i0, i2, i3, i1); // 面(0,2,3)，对面顶点1
            AddFace(faceCount, faceInfo, i1, i2, i3, i0); // 面(1,2,3)，对面顶点0
        }

        var surfTriList = new List<int>();
        foreach (var kvp in faceCount)
        {
            if (kvp.Value == 1)
            {
                var info = faceInfo[kvp.Key];
                int a = info.Item1;
                int b = info.Item2;
                int c = info.Item3;
                int opposite = info.Item4;

                // 确保法线朝外：法线应指向远离对面顶点的方向
                Vector3 pa = vertices[a];
                Vector3 pb = vertices[b];
                Vector3 pc = vertices[c];
                Vector3 po = vertices[opposite];

                Vector3 faceNormal = Vector3.Cross(pb - pa, pc - pa);
                Vector3 toOpposite = po - pa;

                // 如果法线指向对面顶点（朝内），则翻转绕序
                if (Vector3.Dot(faceNormal, toOpposite) > 0f)
                {
                    surfTriList.Add(a);
                    surfTriList.Add(c);
                    surfTriList.Add(b);
                }
                else
                {
                    surfTriList.Add(a);
                    surfTriList.Add(b);
                    surfTriList.Add(c);
                }
            }
        }

        surfaceTriangles = surfTriList.ToArray();
    }

    static void AddEdge(HashSet<(int, int)> set, int a, int b)
    {
        var key = (Mathf.Min(a, b), Mathf.Max(a, b));
        set.Add(key);
    }

    static void AddFace(Dictionary<(int, int, int), int> countDict,
        Dictionary<(int, int, int), (int, int, int, int)> infoDict,
        int a, int b, int c, int opposite)
    {
        // 排序后作为key（用于匹配相邻四面体的共享面）
        int[] sorted = { a, b, c };
        System.Array.Sort(sorted);
        var key = (sorted[0], sorted[1], sorted[2]);

        if (countDict.ContainsKey(key))
        {
            countDict[key]++;
        }
        else
        {
            countDict[key] = 1;
            infoDict[key] = (a, b, c, opposite);
        }
    }

    /// <summary>
    /// 计算四面体有符号体积
    /// V = dot(p1-p0, cross(p2-p0, p3-p0)) / 6
    /// </summary>
    static float ComputeTetVolume(Vector3 p0, Vector3 p1, Vector3 p2, Vector3 p3)
    {
        Vector3 d1 = p1 - p0;
        Vector3 d2 = p2 - p0;
        Vector3 d3 = p3 - p0;
        return Vector3.Dot(d1, Vector3.Cross(d2, d3)) / 6f;
    }

    /// <summary>
    /// 生成球面映射UV
    /// </summary>
    Vector2[] GenerateSphereUVs(Vector3[] verts)
    {
        var uvs = new Vector2[verts.Length];
        Vector3 center = Vector3.zero;
        for (int i = 0; i < verts.Length; i++)
            center += verts[i];
        center /= verts.Length;

        for (int i = 0; i < verts.Length; i++)
        {
            Vector3 dir = (verts[i] - center).normalized;
            float u = 0.5f + Mathf.Atan2(dir.z, dir.x) / (2f * Mathf.PI);
            float v = 0.5f + Mathf.Asin(Mathf.Clamp(dir.y, -1f, 1f)) / Mathf.PI;
            uvs[i] = new Vector2(u, v);
        }
        return uvs;
    }

    /// <summary>
    /// 从顶点、三角形、UV数据创建Mesh
    /// </summary>
    Mesh CreateSurfaceMeshFromData(Vector3[] verts, int[] triangles, Vector2[] uvs)
    {
        Mesh mesh = new Mesh();
        mesh.name = "SoftBody_DOTS_Runtime";

        // 如果顶点数超过65535，使用32位索引
        if (verts.Length > 65535)
            mesh.indexFormat = UnityEngine.Rendering.IndexFormat.UInt32;

        mesh.SetVertices(new List<Vector3>(verts));
        mesh.SetTriangles(triangles, 0);
        if (uvs != null)
            mesh.uv = uvs;

        mesh.RecalculateNormals();
        mesh.RecalculateTangents();
        mesh.RecalculateBounds();
        return mesh;
    }
}
