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
        // 创建Mesh
        clothMesh = CreateClothMesh();

        // 设置渲染
        var meshFilter = gameObject.AddComponent<MeshFilter>();
        var meshRenderer = gameObject.AddComponent<MeshRenderer>();
        meshFilter.mesh = clothMesh;
        meshRenderer.material = clothMaterial;

        var meshVertices = clothMesh.vertices;
        var meshTriangles = clothMesh.triangles;
        int numParticles = meshVertices.Length;

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
        var archetype = entityManager.CreateArchetype(
            typeof(ClothTag),
            typeof(MeshUpdateTag),
            typeof(ClothSolverConfig),
            typeof(ParticlePosition),
            typeof(ParticlePrevPosition),
            typeof(ParticleVelocity),
            typeof(ParticleInvMass),
            typeof(ClothEdge),
            typeof(ClothDistanceLambda),
            typeof(TriangleIndex)
        );

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

        // 三角形索引
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
