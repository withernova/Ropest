using System.Collections.Generic;
using Unity.Collections;
using Unity.Entities;
using Unity.Mathematics;
using UnityEngine;

/// <summary>
/// 运行时Cloth初始化器 - 不依赖SubScene，可以在运行时动态创建Cloth Entity
/// 使用EntityArchetype一次性创建所有组件，避免多次结构性变更导致Buffer失效
/// </summary>
public class ClothRuntimeSpawner : MonoBehaviour
{
    [Header("环境参数")]
    public Vector3 gravity = new(0, -9.8f, 0);

    [Header("布料参数")]
    [Range(0.1f, 100f)] public float length = 5;
    [Range(0.1f, 100f)] public float width = 5;
    [Range(1, 100)] public int segments = 10;
    [Range(1, 100)] public int subdivision = 11;
    public Vector3 meshOrigin;
    public Material clothMaterial;

    [Header("求解器参数")]
    [Range(1, 20)] public int numSubSteps = 7;
    [Range(0f, 1f)] public float distanceStiffness = 0f;
    [Range(0f, 1f)] public float damping = 0.03f;
    [Range(0.01f, 0.5f)] public float collisionRadius = 0.05f;
    [Range(0f, 1f)] public float friction = 0.4f;

    [Header("渲染网格细分")]
    [Tooltip("渲染网格细分迭代次数（0=不细分直接用模拟网格，1=4倍面数，2=16倍面数）")]
    [Range(0, 3)] public int renderSubdivisionIterations = 1;

    [Header("固定点")]
    [Tooltip("固定的顶点索引列表（invMass设为0）")]
    public List<int> fixedVertices = new List<int>();

    private Entity clothEntity;
    private EntityManager entityManager;
    private Mesh clothMesh;
    private bool entityCreated = false;

    void Start()
    {
        entityManager = World.DefaultGameObjectInjectionWorld.EntityManager;
        SpawnCloth();
    }

    void SpawnCloth()
    {
        // 创建模拟用Mesh（低面数）
        var simMesh = CreateClothMesh();
        var meshVertices = simMesh.vertices;
        var meshTriangles = simMesh.triangles;
        var meshUVs = simMesh.uv;
        int numParticles = meshVertices.Length;

        // 生成渲染用Mesh（支持细分高面数）
        bool useSubdivision = renderSubdivisionIterations > 0;
        Vector3[] renderVertices;
        int[] renderTriangles;
        Vector2[] renderUVs;
        SubdivisionMeshHelper.BindingData[] bindingData = null;

        if (useSubdivision)
        {
            // 细分的同时精确追踪每个顶点的重心坐标绑定（无搜索误差）
            SubdivisionMeshHelper.SubdivideWithBindings(
                meshVertices, meshTriangles, meshUVs, renderSubdivisionIterations,
                out renderVertices, out renderTriangles, out renderUVs, out bindingData);

            Debug.Log($"[Cloth] 细分渲染网格: 模拟顶点={numParticles}, 渲染顶点={renderVertices.Length}, " +
                      $"模拟三角形={meshTriangles.Length / 3}, 渲染三角形={renderTriangles.Length / 3}");
        }
        else
        {
            renderVertices = meshVertices;
            renderTriangles = meshTriangles;
            renderUVs = meshUVs;
        }

        // 创建渲染Mesh
        clothMesh = new Mesh();
        clothMesh.name = "Cloth_DOTS_Runtime";
        clothMesh.MarkDynamic(); // 每帧更新顶点，标记为Dynamic VBO
        if (renderVertices.Length > 65535)
            clothMesh.indexFormat = UnityEngine.Rendering.IndexFormat.UInt32;
        clothMesh.SetVertices(renderVertices);
        clothMesh.SetTriangles(renderTriangles, 0);
        if (renderUVs != null) clothMesh.uv = renderUVs;
        clothMesh.RecalculateNormals();
        clothMesh.RecalculateTangents();
        clothMesh.RecalculateBounds();

        // 设置渲染
        var meshFilter = gameObject.AddComponent<MeshFilter>();
        var meshRenderer = gameObject.AddComponent<MeshRenderer>();
        meshFilter.mesh = clothMesh;
        meshRenderer.material = clothMaterial;

        // 计算逆质量
        var invMasses = new float[numParticles];
        for (int i = 0; i < numParticles; i++) invMasses[i] = 1f;

        for (int i = 0; i < meshTriangles.Length; i += 3)
        {
            int i0 = meshTriangles[i];
            int i1 = meshTriangles[i + 1];
            int i2 = meshTriangles[i + 2];

            Vector3 p0 = meshVertices[i0];
            Vector3 p1 = meshVertices[i1];
            Vector3 p2 = meshVertices[i2];

            float a = (p0 - p1).magnitude;
            float b = (p2 - p1).magnitude;
            float c = (p0 - p2).magnitude;
            float p = (a + b + c) / 2f;
            float s = Mathf.Sqrt(p * (p - a) * (p - b) * (p - c));

            invMasses[i0] += s / 3f;
            invMasses[i1] += s / 3f;
            invMasses[i2] += s / 3f;
        }

        for (int i = 0; i < numParticles; i++)
        {
            invMasses[i] = 1f / invMasses[i];
        }

        foreach (int idx in fixedVertices)
        {
            if (idx >= 0 && idx < numParticles)
                invMasses[idx] = 0f;
        }

        // 提取边
        var edgeDict = new Dictionary<(int, int), float>();
        for (int i = 0; i < meshTriangles.Length; i += 3)
        {
            int[] tri = { meshTriangles[i], meshTriangles[i + 1], meshTriangles[i + 2] };
            System.Array.Sort(tri);
            AddEdge(edgeDict, tri[0], tri[1], meshVertices);
            AddEdge(edgeDict, tri[0], tri[2], meshVertices);
            AddEdge(edgeDict, tri[1], tri[2], meshVertices);
        }

        // === 使用Archetype一次性创建Entity ===
        var componentTypes = new List<ComponentType>
        {
            typeof(ClothTag),
            typeof(MeshUpdateTag),
            typeof(ClothSolverConfig),
            typeof(ParticlePosition),
            typeof(ParticlePrevPosition),
            typeof(ParticleVelocity),
            typeof(ParticleInvMass),
            typeof(ClothEdge),
            typeof(ClothDistanceLambda),
            typeof(TriangleIndex),
            typeof(RenderMeshConfig)
        };
        if (useSubdivision)
        {
            componentTypes.Add(typeof(RenderVertexBinding));
        }
        var archetype = entityManager.CreateArchetype(componentTypes.ToArray());

        clothEntity = entityManager.CreateEntity(archetype);

        // 设置配置
        entityManager.SetComponentData(clothEntity, new ClothSolverConfig
        {
            NumParticles = numParticles,
            Subdivision = subdivision,
            NumSubSteps = numSubSteps,
            Gravity = gravity,
            DistanceStiffness = distanceStiffness,
            Damping = damping,
            CollisionRadius = collisionRadius,
            Friction = friction
        });

        // === 填充Buffer数据（GetBuffer不触发结构性变更） ===
        var posBuf = entityManager.GetBuffer<ParticlePosition>(clothEntity);
        var prevBuf = entityManager.GetBuffer<ParticlePrevPosition>(clothEntity);
        var velBuf = entityManager.GetBuffer<ParticleVelocity>(clothEntity);
        var massBuf = entityManager.GetBuffer<ParticleInvMass>(clothEntity);

        for (int i = 0; i < numParticles; i++)
        {
            float3 pos = meshVertices[i];
            posBuf.Add(new ParticlePosition { Value = pos });
            prevBuf.Add(new ParticlePrevPosition { Value = pos });
            velBuf.Add(new ParticleVelocity { Value = float3.zero });
            massBuf.Add(new ParticleInvMass { Value = invMasses[i] });
        }

        // 边数据
        var edgeBuf = entityManager.GetBuffer<ClothEdge>(clothEntity);
        var lambdaBuf = entityManager.GetBuffer<ClothDistanceLambda>(clothEntity);

        foreach (var kvp in edgeDict)
        {
            edgeBuf.Add(new ClothEdge
            {
                IndexA = kvp.Key.Item1,
                IndexB = kvp.Key.Item2,
                RestLength = kvp.Value
            });
            lambdaBuf.Add(new ClothDistanceLambda { Value = 0f });
        }

        // 三角形索引（模拟三角形，用于绑定查找）
        var triBuf = entityManager.GetBuffer<TriangleIndex>(clothEntity);
        for (int i = 0; i < meshTriangles.Length; i += 3)
        {
            triBuf.Add(new TriangleIndex
            {
                I0 = meshTriangles[i],
                I1 = meshTriangles[i + 1],
                I2 = meshTriangles[i + 2]
            });
        }

        // 设置渲染网格配置
        entityManager.SetComponentData(clothEntity, new RenderMeshConfig
        {
            NumRenderVertices = renderVertices.Length,
            UseSubdivision = useSubdivision
        });

        // 填充渲染顶点绑定数据
        if (useSubdivision && bindingData != null)
        {
            var bindBuf = entityManager.GetBuffer<RenderVertexBinding>(clothEntity);
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

        // 托管Mesh引用（最后添加，唯一的额外结构性变更）
        entityManager.AddComponentObject(clothEntity, new ManagedMeshReference
        {
            Mesh = clothMesh,
            MeshFilter = GetComponent<MeshFilter>(),
            MeshRenderer = GetComponent<MeshRenderer>()
        });

        entityCreated = true;
    }

    void OnDestroy()
    {
        if (!entityCreated) return;
        if (World.DefaultGameObjectInjectionWorld == null || !World.DefaultGameObjectInjectionWorld.IsCreated) return;
        if (entityManager.Exists(clothEntity))
        {
            entityManager.DestroyEntity(clothEntity);
        }
    }

    static void AddEdge(Dictionary<(int, int), float> dict, int a, int b, Vector3[] verts)
    {
        var key = (Mathf.Min(a, b), Mathf.Max(a, b));
        if (!dict.ContainsKey(key))
        {
            dict[key] = (verts[a] - verts[b]).magnitude;
        }
    }

    Mesh CreateClothMesh()
    {
        Mesh mesh = new Mesh();
        mesh.name = "Cloth_DOTS_Runtime";
        mesh.MarkDynamic(); // 每帧更新顶点，标记为Dynamic VBO
        List<Vector3> vertices = new List<Vector3>();
        List<int> triangles = new List<int>();
        List<Vector2> uvs = new List<Vector2>();

        for (int i = 0; i < segments + 1; i++)
        {
            Vector3 point = meshOrigin + (length / segments) * new Vector3(i, 0, 0);
            for (int j = 0; j < subdivision + 1; j++)
            {
                Vector3 point1 = point + (width / subdivision) * new Vector3(0, 0, j);
                vertices.Add(point1);
            }
        }

        for (int i = 0; i < segments; i++)
        {
            for (int j = 0; j < subdivision; j++)
            {
                int nextJ = (j + 1) % (subdivision + 1);
                triangles.Add(i * (subdivision + 1) + j);
                triangles.Add(i * (subdivision + 1) + nextJ);
                triangles.Add((i + 1) * (subdivision + 1) + j);
                triangles.Add((i + 1) * (subdivision + 1) + j);
                triangles.Add(i * (subdivision + 1) + nextJ);
                triangles.Add((i + 1) * (subdivision + 1) + nextJ);
            }
        }

        for (int i = 0; i <= segments; i++)
        {
            for (int j = 0; j <= subdivision; j++)
            {
                uvs.Add(new Vector2(j, i));
            }
        }

        mesh.SetVertices(vertices);
        mesh.SetTriangles(triangles, 0);
        mesh.uv = uvs.ToArray();
        mesh.RecalculateNormals();
        mesh.RecalculateTangents();
        return mesh;
    }
}
