using System.Collections.Generic;
using System.Linq;
using Unity.Collections;
using Unity.Entities;
using Unity.Mathematics;
using UnityEngine;

/// <summary>
/// 运行时Rope初始化器 - 不依赖SubScene，可以在运行时动态创建Rope Entity
/// 使用EntityArchetype一次性创建所有组件，避免多次结构性变更导致Buffer失效
/// </summary>
public class RopeRuntimeSpawner : MonoBehaviour
{
    [Header("环境参数")]
    public Vector3 gravity = new(0, -9.8f, 0);

    [Header("绳子参数")]
    [Range(0.1f, 100f)] public float length = 1;
    [Range(1, 100)] public int segments = 16;
    [Range(3, 30)] public int subdivision = 8;
    [Range(0.01f, 10f)] public float thickness = 0.02f;
    public Vector3 meshOrigin;
    public Material ropeMaterial;

    [Header("求解器参数")]
    [Range(1, 20)] public int numSubSteps = 6;
    [Range(0f, 1f)] public float edgeStiffness = 0.0001f;
    [Range(0f, 1f)] public float bendTwistKs = 0.8f;
    public float ghostDistance = 0.1f;
    public float gravityFactor = 1f;
    [Range(0f, 1f)] public float damping = 0.01f;

    private Entity ropeEntity;
    private EntityManager entityManager;
    private Mesh ropeMesh;
    private bool entityCreated = false;

    void Start()
    {
        entityManager = World.DefaultGameObjectInjectionWorld.EntityManager;
        SpawnRope();
    }

    void SpawnRope()
    {
        // 创建Mesh
        ropeMesh = CreateRopeMesh();

        // 设置渲染
        var meshFilter = gameObject.AddComponent<MeshFilter>();
        var meshRenderer = gameObject.AddComponent<MeshRenderer>();
        meshFilter.mesh = ropeMesh;
        meshRenderer.material = ropeMaterial;

        // 计算截面中心点
        int numSections = (ropeMesh.vertexCount - 2) / subdivision;
        var vertices = ropeMesh.vertices;
        var pointPositions = new List<float3>();

        for (int i = 0; i < numSections; i++)
        {
            float3 center = float3.zero;
            for (int j = 0; j < subdivision; j++)
            {
                center += (float3)vertices[i * subdivision + 1 + j];
            }
            center /= subdivision;
            pointPositions.Add(center);
        }

        int numPoints = pointPositions.Count;
        int numGhosts = numPoints - 1;

        // 计算静止长度
        var restLengths = new List<float>();
        for (int i = 0; i < numGhosts; i++)
        {
            restLengths.Add(math.length(pointPositions[i] - pointPositions[i + 1]));
        }

        // 计算Ghost点初始位置
        var ghostPositions = new List<float3>();
        for (int i = 0; i < numGhosts; i++)
        {
            float3 mid = (pointPositions[i] + pointPositions[i + 1]) / 2f;
            mid.y += ghostDistance;
            ghostPositions.Add(mid);
        }

        // 计算初始Darboux向量
        var initDarboux = new List<float3>();
        for (int i = 0; i < numPoints - 2; i++)
        {
            float3x3 a = ComputeMaterialFrame(pointPositions[i], pointPositions[i + 1], ghostPositions[i]);
            float3x3 b = ComputeMaterialFrame(pointPositions[i + 1], pointPositions[i + 2], ghostPositions[i + 1]);
            float3 darboux = ComputeDarbouxVector(a, b, restLengths[0] / 2f);
            initDarboux.Add(darboux);
        }

        // === 使用Archetype一次性创建Entity，避免多次结构性变更 ===
        var archetype = entityManager.CreateArchetype(
            typeof(RopeTag),
            typeof(MeshUpdateTag),
            typeof(RopeSolverConfig),
            typeof(ParticlePosition),
            typeof(ParticlePrevPosition),
            typeof(ParticleVelocity),
            typeof(ParticleInvMass),
            typeof(GhostPosition),
            typeof(GhostPrevPosition),
            typeof(GhostVelocity),
            typeof(GhostInvMass),
            typeof(RopeRestLength),
            typeof(InitDarbouxVector),
            typeof(EdgeLambda0),
            typeof(EdgeLambda1),
            typeof(EdgeLambda2),
            typeof(BendTwistLambda),
            typeof(SectionVertexIndex)
        );

        ropeEntity = entityManager.CreateEntity(archetype);

        // 设置配置组件
        entityManager.SetComponentData(ropeEntity, new RopeSolverConfig
        {
            NumPoints = numPoints,
            NumGhostPoints = numGhosts,
            Segments = segments,
            Subdivision = subdivision,
            NumSubSteps = numSubSteps,
            Radius = thickness,
            GhostDistance = ghostDistance,
            GravityFactor = gravityFactor,
            Gravity = gravity,
            EdgeStiffness = edgeStiffness,
            BendTwistKs = bendTwistKs,
            Damping = damping
        });

        // === 填充Buffer数据（此时不会触发结构性变更，因为Buffer已经存在） ===
        float ctrlMass = 0.3f;
        float invMassValue = 1f;

        var posBuf = entityManager.GetBuffer<ParticlePosition>(ropeEntity);
        var prevBuf = entityManager.GetBuffer<ParticlePrevPosition>(ropeEntity);
        var velBuf = entityManager.GetBuffer<ParticleVelocity>(ropeEntity);
        var massBuf = entityManager.GetBuffer<ParticleInvMass>(ropeEntity);

        for (int i = 0; i < numPoints; i++)
        {
            posBuf.Add(new ParticlePosition { Value = pointPositions[i] });
            prevBuf.Add(new ParticlePrevPosition { Value = pointPositions[i] });
            velBuf.Add(new ParticleVelocity { Value = float3.zero });
            massBuf.Add(new ParticleInvMass { Value = (i == 0) ? ctrlMass : invMassValue });
        }

        // Ghost Buffers
        var ghostPosBuf = entityManager.GetBuffer<GhostPosition>(ropeEntity);
        var ghostPrevBuf = entityManager.GetBuffer<GhostPrevPosition>(ropeEntity);
        var ghostVelBuf = entityManager.GetBuffer<GhostVelocity>(ropeEntity);
        var ghostMassBuf = entityManager.GetBuffer<GhostInvMass>(ropeEntity);

        for (int i = 0; i < numGhosts; i++)
        {
            ghostPosBuf.Add(new GhostPosition { Value = ghostPositions[i] });
            ghostPrevBuf.Add(new GhostPrevPosition { Value = ghostPositions[i] });
            ghostVelBuf.Add(new GhostVelocity { Value = float3.zero });
            ghostMassBuf.Add(new GhostInvMass { Value = invMassValue });
        }

        // 静止长度
        var lenBuf = entityManager.GetBuffer<RopeRestLength>(ropeEntity);
        for (int i = 0; i < numGhosts; i++)
        {
            lenBuf.Add(new RopeRestLength { Value = restLengths[i] });
        }

        // 初始Darboux
        var darbouxBuf = entityManager.GetBuffer<InitDarbouxVector>(ropeEntity);
        for (int i = 0; i < initDarboux.Count; i++)
        {
            darbouxBuf.Add(new InitDarbouxVector { Value = initDarboux[i] });
        }

        // Lambda buffers
        var el0 = entityManager.GetBuffer<EdgeLambda0>(ropeEntity);
        var el1 = entityManager.GetBuffer<EdgeLambda1>(ropeEntity);
        var el2 = entityManager.GetBuffer<EdgeLambda2>(ropeEntity);
        for (int i = 0; i < numPoints; i++)
        {
            el0.Add(new EdgeLambda0 { Value = 0f });
            el1.Add(new EdgeLambda1 { Value = 0f });
            el2.Add(new EdgeLambda2 { Value = 0f });
        }

        var btLambda = entityManager.GetBuffer<BendTwistLambda>(ropeEntity);
        for (int i = 0; i < math.max(0, numPoints - 2); i++)
        {
            btLambda.Add(new BendTwistLambda { Value = float3.zero });
        }

        // 截面索引
        var sectionBuf = entityManager.GetBuffer<SectionVertexIndex>(ropeEntity);
        for (int i = 0; i < numSections; i++)
        {
            for (int j = 0; j < subdivision; j++)
            {
                sectionBuf.Add(new SectionVertexIndex { Value = i * subdivision + 1 + j });
            }
        }

        // 托管Mesh引用（这是唯一的结构性变更，但在所有Buffer填充完毕后执行）
        entityManager.AddComponentObject(ropeEntity, new ManagedMeshReference
        {
            Mesh = ropeMesh,
            MeshFilter = GetComponent<MeshFilter>(),
            MeshRenderer = GetComponent<MeshRenderer>()
        });

        entityCreated = true;
    }

    void OnDestroy()
    {
        if (!entityCreated) return;
        if (World.DefaultGameObjectInjectionWorld == null || !World.DefaultGameObjectInjectionWorld.IsCreated) return;
        if (entityManager.Exists(ropeEntity))
        {
            entityManager.DestroyEntity(ropeEntity);
        }
    }

    Mesh CreateRopeMesh()
    {
        Mesh mesh = new Mesh();
        mesh.name = "Rope_DOTS_Runtime";
        List<Vector3> vertices = new List<Vector3>();
        List<int> triangles = new List<int>();
        List<Vector2> uvs = new List<Vector2>();

        for (int i = 0; i < segments + 1; i++)
        {
            Vector3 point = meshOrigin + (length / segments) * new Vector3(i, 0, 0);
            float angleIncrement = 2 * Mathf.PI / subdivision;
            if (i == 0) vertices.Add(point);
            for (int j = 0; j < subdivision; j++)
            {
                float angle = j * angleIncrement;
                float x = point.x;
                float y = point.y + Mathf.Cos(angle) * thickness;
                float z = point.z + Mathf.Sin(angle) * thickness;
                vertices.Add(new Vector3(x, y, z));
            }
            if (i == segments) vertices.Add(point);
        }

        for (int i = 0; i < segments; i++)
        {
            if (i == 0)
            {
                for (int j = 0; j < subdivision; j++)
                {
                    int nextJ = (j + 1) % subdivision;
                    triangles.Add(0);
                    triangles.Add(1 + nextJ);
                    triangles.Add(1 + j);
                }
            }
            for (int j = 0; j < subdivision; j++)
            {
                int nextJ = (j + 1) % subdivision;
                triangles.Add(i * subdivision + 1 + j);
                triangles.Add(i * subdivision + 1 + nextJ);
                triangles.Add((i + 1) * subdivision + 1 + j);
                triangles.Add((i + 1) * subdivision + 1 + j);
                triangles.Add(i * subdivision + 1 + nextJ);
                triangles.Add((i + 1) * subdivision + 1 + nextJ);
            }
        }
        for (int j = 0; j < subdivision; j++)
        {
            int nextJ = (j + 1) % subdivision;
            triangles.Add(segments * subdivision + 1 + j);
            triangles.Add(segments * subdivision + 1 + nextJ);
            triangles.Add((segments + 1) * subdivision + 1);
        }

        uvs.Add(new Vector2(0, 0));
        for (int i = 0; i <= segments; i++)
        {
            for (int j = 0; j < subdivision; j++)
            {
                float u = (j <= subdivision / 2) ? j / subdivision * 2 : (subdivision - j) / subdivision * 2;
                uvs.Add(new Vector2(u, i));
            }
        }
        uvs.Add(new Vector2(0, 0));

        mesh.SetVertices(vertices);
        mesh.SetTriangles(triangles, 0);
        mesh.uv = uvs.ToArray();
        mesh.RecalculateNormals();
        mesh.RecalculateTangents();
        return mesh;
    }

    static float3x3 ComputeMaterialFrame(float3 p1, float3 p2, float3 pg)
    {
        float3 d3 = math.normalize(p2 - p1);
        float3 d2 = math.normalize(math.cross(d3, pg - p1));
        float3 d1 = math.cross(d2, d3);
        return new float3x3(d1, d2, d3);
    }

    static float3 ComputeDarbouxVector(float3x3 dA, float3x3 dB, float halfLength)
    {
        float factor = 1.0f + math.dot(dA.c0, dB.c0) + math.dot(dA.c1, dB.c1) + math.dot(dA.c2, dB.c2);
        factor = 2.0f / (halfLength * factor);

        float3 darboux;
        darboux.x = math.dot(dA.c2, dB.c1) - math.dot(dA.c1, dB.c2);
        darboux.y = math.dot(dA.c0, dB.c2) - math.dot(dA.c2, dB.c0);
        darboux.z = math.dot(dA.c1, dB.c0) - math.dot(dA.c0, dB.c1);

        return darboux * factor;
    }
}
