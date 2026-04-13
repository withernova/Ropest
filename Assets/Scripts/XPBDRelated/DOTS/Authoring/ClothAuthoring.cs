using System.Collections.Generic;
using System.Linq;
using Unity.Collections;
using Unity.Entities;
using Unity.Mathematics;
using UnityEngine;

/// <summary>
/// Cloth Authoring组件 - 挂载到GameObject上，Baker会将其转换为ECS Entity
/// 替代原有的Cloth MonoBehaviour
/// </summary>
public class ClothAuthoring : MonoBehaviour
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

    [Header("固定点")]
    [Tooltip("固定的顶点索引列表（invMass设为0）")]
    public List<int> fixedVertices = new List<int>();

    class Baker : Baker<ClothAuthoring>
    {
        public override void Bake(ClothAuthoring authoring)
        {
            var entity = GetEntity(TransformUsageFlags.Dynamic);

            // 创建Mesh
            Mesh clothMesh = CreateClothMesh(authoring);
            var meshVertices = clothMesh.vertices;
            var meshTriangles = clothMesh.triangles;
            int numParticles = meshVertices.Length;

            // 计算逆质量
            var invMasses = new float[numParticles];
            for (int i = 0; i < numParticles; i++) invMasses[i] = 1f;

            // 计算三角形面积并分配质量
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

            // 固定点
            foreach (int idx in authoring.fixedVertices)
            {
                if (idx >= 0 && idx < numParticles)
                    invMasses[idx] = 0f;
            }

            // 提取边（距离约束）
            var edgeDict = new Dictionary<(int, int), float>();
            for (int i = 0; i < meshTriangles.Length; i += 3)
            {
                int[] tri = { meshTriangles[i], meshTriangles[i + 1], meshTriangles[i + 2] };
                System.Array.Sort(tri);

                AddEdge(edgeDict, tri[0], tri[1], meshVertices);
                AddEdge(edgeDict, tri[0], tri[2], meshVertices);
                AddEdge(edgeDict, tri[1], tri[2], meshVertices);
            }

            // === 添加组件 ===
            AddComponent(entity, new ClothTag());
            AddComponent(entity, new MeshUpdateTag());
            AddComponent(entity, new ClothSolverConfig
            {
                NumParticles = numParticles,
                Subdivision = authoring.subdivision,
                NumSubSteps = authoring.numSubSteps,
                Gravity = authoring.gravity,
                DistanceStiffness = authoring.distanceStiffness
            });

            // === 添加Buffers ===
            var posBuf = AddBuffer<ParticlePosition>(entity);
            var prevBuf = AddBuffer<ParticlePrevPosition>(entity);
            var velBuf = AddBuffer<ParticleVelocity>(entity);
            var massBuf = AddBuffer<ParticleInvMass>(entity);

            for (int i = 0; i < numParticles; i++)
            {
                float3 pos = meshVertices[i];
                posBuf.Add(new ParticlePosition { Value = pos });
                prevBuf.Add(new ParticlePrevPosition { Value = pos });
                velBuf.Add(new ParticleVelocity { Value = float3.zero });
                massBuf.Add(new ParticleInvMass { Value = invMasses[i] });
            }

            // 边数据
            var edgeBuf = AddBuffer<ClothEdge>(entity);
            var lambdaBuf = AddBuffer<ClothDistanceLambda>(entity);

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
            var triBuf = AddBuffer<TriangleIndex>(entity);
            for (int i = 0; i < meshTriangles.Length; i += 3)
            {
                triBuf.Add(new TriangleIndex
                {
                    I0 = meshTriangles[i],
                    I1 = meshTriangles[i + 1],
                    I2 = meshTriangles[i + 2]
                });
            }

            // 托管Mesh引用
            AddComponentObject(entity, new ManagedMeshReference
            {
                Mesh = clothMesh
            });
        }

        static void AddEdge(Dictionary<(int, int), float> dict, int a, int b, Vector3[] verts)
        {
            var key = (Mathf.Min(a, b), Mathf.Max(a, b));
            if (!dict.ContainsKey(key))
            {
                dict[key] = (verts[a] - verts[b]).magnitude;
            }
        }

        static Mesh CreateClothMesh(ClothAuthoring authoring)
        {
            float len = authoring.length;
            float wid = authoring.width;
            int seg = authoring.segments;
            int sub = authoring.subdivision;
            Vector3 initPos = authoring.meshOrigin;

            Mesh mesh = new Mesh();
            mesh.name = "Cloth_DOTS";
            List<Vector3> vertices = new List<Vector3>();
            List<int> triangles = new List<int>();
            List<Vector2> uvs = new List<Vector2>();

            // 顶点
            for (int i = 0; i < seg + 1; i++)
            {
                Vector3 point = initPos + (len / seg) * new Vector3(i, 0, 0);
                for (int j = 0; j < sub + 1; j++)
                {
                    Vector3 point1 = point + (wid / sub) * new Vector3(0, 0, j);
                    vertices.Add(point1);
                }
            }

            // 索引
            for (int i = 0; i < seg; i++)
            {
                for (int j = 0; j < sub; j++)
                {
                    int nextJ = (j + 1) % (sub + 1);
                    triangles.Add(i * (sub + 1) + j);
                    triangles.Add(i * (sub + 1) + nextJ);
                    triangles.Add((i + 1) * (sub + 1) + j);

                    triangles.Add((i + 1) * (sub + 1) + j);
                    triangles.Add(i * (sub + 1) + nextJ);
                    triangles.Add((i + 1) * (sub + 1) + nextJ);
                }
            }

            // uv
            for (int i = 0; i <= seg; i++)
            {
                for (int j = 0; j <= sub; j++)
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
}
