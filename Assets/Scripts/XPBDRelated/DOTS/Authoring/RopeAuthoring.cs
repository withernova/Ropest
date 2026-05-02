using System.Collections.Generic;
using System.Linq;
using Unity.Collections;
using Unity.Entities;
using Unity.Mathematics;
using UnityEngine;

public class RopeAuthoring : MonoBehaviour
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
    // XPBD compliance (α)，值越小越硬。Edge 原版默认 0.0001（几乎刚性），BendTwist 原版默认 0.6（较软，用来稳定弯扭）。
    [Range(0f, 1f)] public float edgeStiffness = 0.0001f;
    [Range(0f, 10f)] public float bendTwistStiffness = 0.6f;
    // BendTwist 的刚度缩放系数 Ks（Position-based Elastic Rods 论文公式中的 Ks），原版默认 0.8。
    [Range(0f, 1f)] public float bendTwistKs = 0.8f;
    public float ghostDistance = 0.1f;
    public float gravityFactor = 1f;
    [Range(0f, 1f)] public float friction = 0.3f;

    class Baker : Baker<RopeAuthoring>
    {
        public override void Bake(RopeAuthoring authoring)
        {
            var entity = GetEntity(TransformUsageFlags.Dynamic);

            // 创建Mesh
            Mesh ropeMesh = CreateRopeMesh(authoring);

            // 计算截面中心点（pointPos）
            int numSections = (ropeMesh.vertexCount - 2) / authoring.subdivision;
            var vertices = ropeMesh.vertices;
            var pointPositions = new List<float3>();

            for (int i = 0; i < numSections; i++)
            {
                float3 center = float3.zero;
                for (int j = 0; j < authoring.subdivision; j++)
                {
                    center += (float3)vertices[i * authoring.subdivision + 1 + j];
                }
                center /= authoring.subdivision;
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
                mid.y += authoring.ghostDistance;
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

            // === 添加组件 ===
            AddComponent(entity, new RopeTag());
            AddComponent(entity, new MeshUpdateTag());
            AddComponent(entity, new RopeSolverConfig
            {
                NumPoints = numPoints,
                NumGhostPoints = numGhosts,
                Segments = authoring.segments,
                Subdivision = authoring.subdivision,
                NumSubSteps = authoring.numSubSteps,
                Radius = authoring.thickness,
                GhostDistance = authoring.ghostDistance,
                GravityFactor = authoring.gravityFactor,
                Gravity = authoring.gravity,
                EdgeStiffness = authoring.edgeStiffness,
                BendTwistStiffness = authoring.bendTwistStiffness,
                BendTwistKs = authoring.bendTwistKs,
                Friction = authoring.friction
            });

            // === 添加Buffers ===
            // 粒子位置
            var posBuf = AddBuffer<ParticlePosition>(entity);
            var prevBuf = AddBuffer<ParticlePrevPosition>(entity);
            var velBuf = AddBuffer<ParticleVelocity>(entity);
            var massBuf = AddBuffer<ParticleInvMass>(entity);

            float invMassValue = 1f;
            float ctrlMass = 0.3f;

            for (int i = 0; i < numPoints; i++)
            {
                posBuf.Add(new ParticlePosition { Value = pointPositions[i] });
                prevBuf.Add(new ParticlePrevPosition { Value = pointPositions[i] });
                velBuf.Add(new ParticleVelocity { Value = float3.zero });
                massBuf.Add(new ParticleInvMass { Value = (i == 0) ? ctrlMass : invMassValue });
            }

            // Ghost点
            var ghostPosBuf = AddBuffer<GhostPosition>(entity);
            var ghostPrevBuf = AddBuffer<GhostPrevPosition>(entity);
            var ghostVelBuf = AddBuffer<GhostVelocity>(entity);
            var ghostMassBuf = AddBuffer<GhostInvMass>(entity);

            for (int i = 0; i < numGhosts; i++)
            {
                ghostPosBuf.Add(new GhostPosition { Value = ghostPositions[i] });
                ghostPrevBuf.Add(new GhostPrevPosition { Value = ghostPositions[i] });
                ghostVelBuf.Add(new GhostVelocity { Value = float3.zero });
                ghostMassBuf.Add(new GhostInvMass { Value = invMassValue });
            }

            // 静止长度
            var lenBuf = AddBuffer<RopeRestLength>(entity);
            for (int i = 0; i < numGhosts; i++)
            {
                lenBuf.Add(new RopeRestLength { Value = restLengths[i] });
            }

            // 初始Darboux
            var darbouxBuf = AddBuffer<InitDarbouxVector>(entity);
            for (int i = 0; i < initDarboux.Count; i++)
            {
                darbouxBuf.Add(new InitDarbouxVector { Value = initDarboux[i] });
            }

            // Lambda buffers
            var el0 = AddBuffer<EdgeLambda0>(entity);
            var el1 = AddBuffer<EdgeLambda1>(entity);
            var el2 = AddBuffer<EdgeLambda2>(entity);
            for (int i = 0; i < numPoints; i++)
            {
                el0.Add(new EdgeLambda0 { Value = 0f });
                el1.Add(new EdgeLambda1 { Value = 0f });
                el2.Add(new EdgeLambda2 { Value = 0f });
            }

            var btLambda = AddBuffer<BendTwistLambda>(entity);
            for (int i = 0; i < math.max(0, numPoints - 2); i++)
            {
                btLambda.Add(new BendTwistLambda { Value = float3.zero });
            }

            // 截面顶点索引（用于渲染映射）
            var sectionBuf = AddBuffer<SectionVertexIndex>(entity);
            for (int i = 0; i < numSections; i++)
            {
                for (int j = 0; j < authoring.subdivision; j++)
                {
                    sectionBuf.Add(new SectionVertexIndex { Value = i * authoring.subdivision + 1 + j });
                }
            }

            // 托管Mesh引用（通过ManagedMeshReference）
            AddComponentObject(entity, new ManagedMeshReference
            {
                Mesh = ropeMesh
            });
        }

        static Mesh CreateRopeMesh(RopeAuthoring authoring)
        {
            float len = authoring.length;
            int seg = authoring.segments;
            int sub = authoring.subdivision;
            float radius = authoring.thickness;
            Vector3 initPos = authoring.meshOrigin;

            Mesh mesh = new Mesh();
            mesh.name = "Rope_DOTS";
            mesh.MarkDynamic();
            List<Vector3> vertices = new List<Vector3>();
            List<int> triangles = new List<int>();
            List<Vector2> uvs = new List<Vector2>();

            // 顶点
            for (int i = 0; i < seg + 1; i++)
            {
                Vector3 point = initPos + (len / seg) * new Vector3(i, 0, 0);
                float angleIncrement = 2 * Mathf.PI / sub;
                if (i == 0)
                    vertices.Add(point);
                for (int j = 0; j < sub; j++)
                {
                    float angle = j * angleIncrement;
                    float x = point.x;
                    float y = point.y + Mathf.Cos(angle) * radius;
                    float z = point.z + Mathf.Sin(angle) * radius;
                    vertices.Add(new Vector3(x, y, z));
                }
                if (i == seg)
                    vertices.Add(point);
            }

            // 索引
            for (int i = 0; i < seg; i++)
            {
                if (i == 0)
                {
                    for (int j = 0; j < sub; j++)
                    {
                        int nextJ = (j + 1) % sub;
                        triangles.Add(0);
                        triangles.Add(1 + nextJ);
                        triangles.Add(1 + j);
                    }
                }
                for (int j = 0; j < sub; j++)
                {
                    int nextJ = (j + 1) % sub;
                    triangles.Add(i * sub + 1 + j);
                    triangles.Add(i * sub + 1 + nextJ);
                    triangles.Add((i + 1) * sub + 1 + j);

                    triangles.Add((i + 1) * sub + 1 + j);
                    triangles.Add(i * sub + 1 + nextJ);
                    triangles.Add((i + 1) * sub + 1 + nextJ);
                }
            }
            for (int j = 0; j < sub; j++)
            {
                int nextJ = (j + 1) % sub;
                triangles.Add(seg * sub + 1 + j);
                triangles.Add(seg * sub + 1 + nextJ);
                triangles.Add((seg + 1) * sub + 1);
            }

            // uv
            uvs.Add(new Vector2(0, 0));
            for (int i = 0; i <= seg; i++)
            {
                for (int j = 0; j < sub; j++)
                {
                    float u, v;
                    if (j <= sub / 2) u = j / sub * 2;
                    else u = (sub - j) / sub * 2;
                    v = i;
                    uvs.Add(new Vector2(u, v));
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
}
