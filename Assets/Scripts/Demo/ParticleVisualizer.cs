using System.Collections.Generic;
using Unity.Collections;
using Unity.Entities;
using Unity.Mathematics;
using UnityEngine;

namespace Ropest.Demo
{
    /// <summary>
    /// 粒子可视化阶段渲染器（URP 兼容版）。
    ///
    /// 上一版用 OnRenderObject + GL.LINES/QUADS，那是 BuiltIn 时代的 API；
    /// 在 URP 下 OnRenderObject 不会被调用，所以画不出来。
    ///
    /// 这一版改成"每帧重建两个 Mesh（点用 Quad、线用 Lines），用 Graphics.DrawMesh 提交"：
    ///   - DrawMesh 走的是正常 URP 渲染管线，URP / BuiltIn 都能正确显示；
    ///   - 用 MeshTopology.Lines 让每两个顶点构成一段线（速度向量）；
    ///   - 粒子点用 billboard quad（朝相机），单段几千顶点完全没压力。
    ///
    /// 材质：直接用 URP 的 Unlit 或 BuiltIn 的 Sprites/Default 都行，这里用 Sprites/Default
    ///       因为它对 vertex color 友好且一定存在。
    /// </summary>
    [DisallowMultipleComponent]
    public class ParticleVisualizer : MonoBehaviour
    {
        [Header("绘制开关")]
        public bool DrawParticles = true;
        public bool DrawVelocityLines = true;

        [Header("粒子点样式")]
        [Tooltip("粒子点（billboard 圆点）的世界半径")]
        public float ParticleSize = 0.018f;
        public Color ParticleColor = new Color(0.10f, 0.10f, 0.12f, 1f);

        [Header("速度线样式")]
        [Tooltip("速度向量缩放：实际线长 = velocity * scale")]
        public float VelocityScale = 0.08f;
        [Tooltip("低于该速度的粒子不画线（避免静止粒子堆满屏）")]
        public float MinSpeedToDrawLine = 0.05f;
        public Color VelocityLineColor = new Color(1f, 0.95f, 0.15f, 1f);   // 亮黄
        public Color VelocityHeadColor = new Color(1f, 0.55f, 0.05f, 1f);   // 橙

        // 我们隐藏过的 MeshRenderer 集合（用于 Disable 时还原；
        // 用 HashSet 是为了 RescanAndHideNew 做"已隐藏过"快速判定，避免重复处理同一个 MR）
        readonly HashSet<MeshRenderer> _hiddenRenderers = new HashSet<MeshRenderer>();

        Camera _targetCamera;

        // 动态 Mesh：每帧重建。Mesh 用 UInt32 索引，最多支持几十万顶点。
        Mesh _pointMesh;
        Mesh _lineMesh;

        // 点和线用两套材质：
        //   - 点：Sprites/Default + 圆形 alpha 贴图，让 quad 看起来是圆点（带边缘渐隐，更柔和）
        //   - 线：Sprites/Default 不带贴图（_MainTex = 白），纯顶点色显示
        // 不能共用一个材质：线的顶点没有 UV，会采样到贴图边缘的 alpha=0 区域而消失。
        Material _pointMaterial;
        Material _lineMaterial;
        Texture2D _circleTexture;

        // 复用的 CPU 缓冲，避免每帧 GC
        readonly List<Vector3> _vertBuf = new List<Vector3>(8192);
        readonly List<Color> _colorBuf = new List<Color>(8192);
        readonly List<Vector2> _uvBuf = new List<Vector2>(8192);
        readonly List<int> _indexBuf = new List<int>(8192);

        public void Activate(Camera targetCamera)
        {
            _targetCamera = targetCamera;
            HideSimulationMeshes();
            EnsureResources();
        }

        public void Deactivate()
        {
            RestoreSimulationMeshes();
        }

        void OnDestroy()
        {
            RestoreSimulationMeshes();
            if (_pointMesh != null) Destroy(_pointMesh);
            if (_lineMesh != null) Destroy(_lineMesh);
            if (_pointMaterial != null) Destroy(_pointMaterial);
            if (_lineMaterial != null) Destroy(_lineMaterial);
            if (_circleTexture != null) Destroy(_circleTexture);
        }

        void EnsureResources()
        {
            if (_pointMesh == null)
            {
                _pointMesh = new Mesh { name = "ParticleVis_Points" };
                _pointMesh.indexFormat = UnityEngine.Rendering.IndexFormat.UInt32;
                _pointMesh.MarkDynamic();
            }
            if (_lineMesh == null)
            {
                _lineMesh = new Mesh { name = "ParticleVis_Lines" };
                _lineMesh.indexFormat = UnityEngine.Rendering.IndexFormat.UInt32;
                _lineMesh.MarkDynamic();
            }
            // 圆形贴图：32×32，中心 alpha=1，向边缘平滑过渡到 0
            if (_circleTexture == null)
            {
                const int s = 32;
                _circleTexture = new Texture2D(s, s, TextureFormat.Alpha8, false)
                {
                    name = "ParticleVis_CircleTex",
                    hideFlags = HideFlags.HideAndDontSave,
                    wrapMode = TextureWrapMode.Clamp,
                    filterMode = FilterMode.Bilinear
                };
                var pixels = new Color32[s * s];
                float center = (s - 1) * 0.5f;
                float maxR = center;
                for (int y = 0; y < s; y++)
                {
                    for (int x = 0; x < s; x++)
                    {
                        float dx = x - center;
                        float dy = y - center;
                        float r = Mathf.Sqrt(dx * dx + dy * dy);
                        // 0..0.85 之间 alpha=1，0.85..1.0 平滑过渡到 0
                        float t = Mathf.InverseLerp(maxR * 0.85f, maxR, r);
                        float a = 1f - t;
                        a = Mathf.Clamp01(a);
                        pixels[y * s + x] = new Color32(255, 255, 255, (byte)(a * 255f));
                    }
                }
                _circleTexture.SetPixels32(pixels);
                _circleTexture.Apply(false, true);
            }

            var shader = Shader.Find("Sprites/Default");
            if (shader == null) shader = Shader.Find("Unlit/Color"); // fallback

            if (_pointMaterial == null)
            {
                _pointMaterial = new Material(shader)
                {
                    hideFlags = HideFlags.HideAndDontSave,
                    name = "ParticleVis_PointMat"
                };
                _pointMaterial.mainTexture = _circleTexture;
            }
            if (_lineMaterial == null)
            {
                _lineMaterial = new Material(shader)
                {
                    hideFlags = HideFlags.HideAndDontSave,
                    name = "ParticleVis_LineMat"
                };
                // 不设 mainTexture，Sprites/Default 默认采样到全白，乘以顶点色 = 顶点色本身
            }
        }

        void HideSimulationMeshes()
        {
            // 注意：可视化阶段是"重播"，仿真物体在 Activate 之后才被段的 OnEnter 重新生成，
            // 所以这里第一次扫到的可能是 0 个；真正的隐藏发生在 LateUpdate 里的 RescanAndHideNew。
#if UNITY_2023_1_OR_NEWER
            var allRenderers = Object.FindObjectsByType<MeshRenderer>(FindObjectsSortMode.None);
#else
            var allRenderers = Object.FindObjectsOfType<MeshRenderer>();
#endif
            foreach (var mr in allRenderers)
            {
                TryHideIfSimMesh(mr);
            }
        }

        /// <summary>每帧重扫一次：捕获重播阶段新生成的仿真 MeshRenderer 并隐藏。</summary>
        void RescanAndHideNew()
        {
#if UNITY_2023_1_OR_NEWER
            var allRenderers = Object.FindObjectsByType<MeshRenderer>(FindObjectsSortMode.None);
#else
            var allRenderers = Object.FindObjectsOfType<MeshRenderer>();
#endif
            foreach (var mr in allRenderers)
            {
                if (mr == null) continue;
                if (_hiddenRenderers.Contains(mr)) continue; // 已经隐藏过的不重复处理
                if (!mr.enabled) continue;                    // 还没启用的（如布料延迟启用）下次再说
                TryHideIfSimMesh(mr);
            }
        }

        void TryHideIfSimMesh(MeshRenderer mr)
        {
            if (mr == null) return;
            if (_hiddenRenderers.Contains(mr)) return;
            bool isSimMesh =
                mr.GetComponent<ClothRuntimeSpawner>() != null ||
                mr.GetComponent<SoftBodyRuntimeSpawner>() != null;
            if (!isSimMesh) return;
            mr.enabled = false;
            _hiddenRenderers.Add(mr);
        }

        void RestoreSimulationMeshes()
        {
            foreach (var mr in _hiddenRenderers)
            {
                if (mr != null) mr.enabled = true;
            }
            _hiddenRenderers.Clear();
        }

        void LateUpdate()
        {
            if (_targetCamera == null) return;

            // 持续扫一遍：重播阶段新生成的 Spawner 在自己的 Start() 里才创建 MeshRenderer，
            // 所以必须每帧扫，否则会闪一两帧的实体网格。
            RescanAndHideNew();

            var world = World.DefaultGameObjectInjectionWorld;
            if (world == null || !world.IsCreated) return;
            var em = world.EntityManager;
            // 在读 Buffer 前先 Complete 所有 Job
            em.CompleteAllTrackedJobs();

            if (DrawVelocityLines)
            {
                BuildLineMesh(em);
                if (_lineMesh.vertexCount > 0)
                {
                    Graphics.DrawMesh(_lineMesh, Matrix4x4.identity, _lineMaterial, 0,
                        _targetCamera, 0, null, false, false, false);
                }
            }

            if (DrawParticles)
            {
                BuildPointMesh(em);
                if (_pointMesh.vertexCount > 0)
                {
                    Graphics.DrawMesh(_pointMesh, Matrix4x4.identity, _pointMaterial, 0,
                        _targetCamera, 0, null, false, false, false);
                }
            }
        }

        // ===== 速度线 Mesh：每两个顶点 = 一段线 =====
        void BuildLineMesh(EntityManager em)
        {
            _vertBuf.Clear();
            _colorBuf.Clear();
            _indexBuf.Clear();

            CollectLineVerts<ClothTag>(em);
            CollectLineVerts<SoftBodyTag>(em);

            _lineMesh.Clear();
            if (_vertBuf.Count == 0) return;
            _lineMesh.SetVertices(_vertBuf);
            _lineMesh.SetColors(_colorBuf);
            _lineMesh.SetIndices(_indexBuf, MeshTopology.Lines, 0, false);
            _lineMesh.bounds = new Bounds(Vector3.zero, Vector3.one * 1000f);
        }

        void CollectLineVerts<TTag>(EntityManager em) where TTag : unmanaged, IComponentData
        {
            var query = em.CreateEntityQuery(
                ComponentType.ReadOnly<TTag>(),
                ComponentType.ReadOnly<ParticlePosition>(),
                ComponentType.ReadOnly<ParticleVelocity>());
            var entities = query.ToEntityArray(Allocator.Temp);

            for (int e = 0; e < entities.Length; e++)
            {
                var entity = entities[e];
                var pos = em.GetBuffer<ParticlePosition>(entity, true);
                var vel = em.GetBuffer<ParticleVelocity>(entity, true);
                int n = math.min(pos.Length, vel.Length);

                for (int p = 0; p < n; p++)
                {
                    float3 v = vel[p].Value;
                    float speed = math.length(v);
                    if (speed < MinSpeedToDrawLine) continue;

                    float3 start = pos[p].Value;
                    float3 end = start + v * VelocityScale;

                    int baseIdx = _vertBuf.Count;
                    _vertBuf.Add(new Vector3(start.x, start.y, start.z));
                    _vertBuf.Add(new Vector3(end.x, end.y, end.z));
                    _colorBuf.Add(VelocityLineColor);
                    _colorBuf.Add(VelocityHeadColor);
                    _indexBuf.Add(baseIdx);
                    _indexBuf.Add(baseIdx + 1);
                }
            }
            entities.Dispose();
            query.Dispose();
        }

        // ===== 粒子点 Mesh：每个粒子 = billboard quad（4 顶点 + 2 三角形 + 4 UV）=====
        void BuildPointMesh(EntityManager em)
        {
            _vertBuf.Clear();
            _colorBuf.Clear();
            _uvBuf.Clear();
            _indexBuf.Clear();

            Transform camTrans = _targetCamera.transform;
            Vector3 camRight = camTrans.right;
            Vector3 camUp = camTrans.up;
            float halfSize = ParticleSize;

            CollectPointQuads<ClothTag>(em, camRight, camUp, halfSize);
            CollectPointQuads<SoftBodyTag>(em, camRight, camUp, halfSize);

            _pointMesh.Clear();
            if (_vertBuf.Count == 0) return;
            _pointMesh.SetVertices(_vertBuf);
            _pointMesh.SetColors(_colorBuf);
            _pointMesh.SetUVs(0, _uvBuf);
            _pointMesh.SetIndices(_indexBuf, MeshTopology.Triangles, 0, false);
            _pointMesh.bounds = new Bounds(Vector3.zero, Vector3.one * 1000f);
        }

        void CollectPointQuads<TTag>(EntityManager em, Vector3 camRight, Vector3 camUp, float halfSize)
            where TTag : unmanaged, IComponentData
        {
            var query = em.CreateEntityQuery(
                ComponentType.ReadOnly<TTag>(),
                ComponentType.ReadOnly<ParticlePosition>());
            var entities = query.ToEntityArray(Allocator.Temp);

            for (int e = 0; e < entities.Length; e++)
            {
                var entity = entities[e];
                var pos = em.GetBuffer<ParticlePosition>(entity, true);
                int n = pos.Length;
                for (int p = 0; p < n; p++)
                {
                    float3 c = pos[p].Value;
                    Vector3 center = new Vector3(c.x, c.y, c.z);

                    Vector3 v0 = center + (-camRight - camUp) * halfSize;
                    Vector3 v1 = center + ( camRight - camUp) * halfSize;
                    Vector3 v2 = center + ( camRight + camUp) * halfSize;
                    Vector3 v3 = center + (-camRight + camUp) * halfSize;

                    int baseIdx = _vertBuf.Count;
                    _vertBuf.Add(v0); _colorBuf.Add(ParticleColor); _uvBuf.Add(new Vector2(0f, 0f));
                    _vertBuf.Add(v1); _colorBuf.Add(ParticleColor); _uvBuf.Add(new Vector2(1f, 0f));
                    _vertBuf.Add(v2); _colorBuf.Add(ParticleColor); _uvBuf.Add(new Vector2(1f, 1f));
                    _vertBuf.Add(v3); _colorBuf.Add(ParticleColor); _uvBuf.Add(new Vector2(0f, 1f));

                    _indexBuf.Add(baseIdx + 0);
                    _indexBuf.Add(baseIdx + 1);
                    _indexBuf.Add(baseIdx + 2);
                    _indexBuf.Add(baseIdx + 0);
                    _indexBuf.Add(baseIdx + 2);
                    _indexBuf.Add(baseIdx + 3);
                }
            }
            entities.Dispose();
            query.Dispose();
        }
    }
}
