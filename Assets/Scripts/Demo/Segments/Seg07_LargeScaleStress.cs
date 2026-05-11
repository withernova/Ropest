using UnityEngine;

namespace Ropest.Demo.Segments
{
    /// <summary>
    /// 段 7：大规模性能 + 跨体碰撞合并展示（C 子方案 + 性能段）。
    /// 内容：
    ///   - 1 块 256×256（≈ 65k 粒子）的高分辨率布料，悬挂在四角附近的两条边
    ///   - 上方依次掉落 3 个软体球，砸到布料上弹开
    /// 看点：高分辨率 XPBD 实时性 + 大规模布料与多软体的稳定跨体碰撞。
    /// 镜头：缓慢环绕拉远，让画面"显大"。
    /// </summary>
    public class Seg07_LargeScaleStress : DemoSegmentBase
    {
        public override string Title => "大规模布料 × 多软体（性能演示）";
        public override string Subtitle => "128×128 ≈ 5k 粒子布料 + 3 软体球跨体碰撞";
        public override float Duration => 13f;
        // Seg07 也加可视化重播：与 Duration 一致（基类默认行为，这里不再覆写）。
        // 注意：粒子数极大时可视化可能掉帧，但视觉上"点阵起伏"非常震撼，值得一试。
        public override DemoCameraPreset Camera =>
            DemoCameraPreset.OrbitAndZoom(
                lookAt: new Vector3(0f, 4.0f, 0f),
                distA: 8f, distB: 5f,
                yawA: -35f, yawB: 25f,
                pitch: 40f, fov: 45f);

        // 软体球依次落下的时间表（相对段开始的秒数）
        readonly float[] _ballSpawnTimes = { 1.5f, 4.0f, 6.5f };
        readonly Vector3[] _ballSpawnPositions =
        {
            new Vector3( 0.0f, 16.0f,  0.0f),
            new Vector3(-1.5f, 17.0f,  0.8f),
            new Vector3( 1.4f, 18.0f, -0.6f)
        };
        bool[] _ballSpawned;

        protected override void OnEnter()
        {
            // ===== 高分辨率布料（蹦床式：四角不固定，改为顶端两条边的 4 个角点固定，让中间能托住球）=====
            var clothGo = new GameObject("Seg07_BigCloth");
            clothGo.transform.SetParent(Root, false);
            var sp = clothGo.AddComponent<ClothRuntimeSpawner>();
            sp.length = 6f;
            sp.width = 6f;
            sp.segments = 127;        // 256 列
            sp.subdivision = 127;     // 256 行 → 65536 粒子
            sp.meshOrigin = new Vector3(-3f, 10f, -3f);
            sp.numSubSteps = 6;       // 大规模略降，保证 60fps
            sp.distanceStiffness = 0.0f;
            sp.damping = 0.01f;
            sp.collisionRadius = 0.1f;
            sp.friction = 0.4f;
            sp.useGraphColoring = true;
            sp.enableSelfCollision = false; // 大规模自碰撞太重，关闭
            sp.renderSubdivisionIterations = 0; // 已经够细，不再加细分
            sp.clothMaterial = Director.ClothMat;
            clothGo.AddComponent<ClothRendererEnabler>();

            // 固定四个角，做"被四角拉住的蹦床"
            int subN = sp.subdivision + 1;       // 256
            int segN = sp.segments + 1;          // 256
            int idxCornerA = 0;                                   // (i=0,j=0)
            int idxCornerB = subN - 1;                            // (i=0,j=last)
            int idxCornerC = (segN - 1) * subN;                   // (i=last,j=0)
            int idxCornerD = (segN - 1) * subN + (subN - 1);      // (i=last,j=last)
            sp.fixedVertices = new System.Collections.Generic.List<int>
            {
                idxCornerA, idxCornerB, idxCornerC, idxCornerD
            };

            _ballSpawned = new bool[_ballSpawnTimes.Length];
        }

        protected override void OnUpdate(float dt)
        {
            // 按时间表依次生成 3 个软体球
            for (int i = 0; i < _ballSpawnTimes.Length; i++)
            {
                if (_ballSpawned[i]) continue;
                if (ElapsedTime < _ballSpawnTimes[i]) continue;
                _ballSpawned[i] = true;
                SpawnSoftBall(i, _ballSpawnPositions[i]);
            }
        }

        void SpawnSoftBall(int idx, Vector3 pos)
        {
            var go = new GameObject($"Seg07_SoftBall_{idx}");
            go.transform.SetParent(Root, false);
            var sp = go.AddComponent<SoftBodyRuntimeSpawner>();
            sp.shape = SoftBodyRuntimeSpawner.SoftBodyShape.Sphere;
            float size = 1.0f + idx * 0.1f;
            sp.sizeX = sp.sizeY = sp.sizeZ = size;
            sp.resolution = 7;
            sp.projectToSphereSurface = true;
            sp.meshOrigin = pos - new Vector3(size, size, size) * 0.5f;
            sp.numSubSteps = 10;
            sp.distanceStiffness = 0.0005f;
            sp.volumeStiffness = 0.0000005f;
            sp.damping = 0.02f;
            sp.collisionRadius = 0.1f;
            sp.friction = 0.4f;
            sp.useGraphColoring = true;
            sp.renderSubdivisionIterations = 1;

            // 每个球颜色不同
            var mat = new Material(Director.SoftBodyMat.shader);
            Color[] palette =
            {
                new Color(0.30f, 0.65f, 0.95f, 1f),
                new Color(0.95f, 0.55f, 0.30f, 1f),
                new Color(0.55f, 0.85f, 0.40f, 1f)
            };
            mat.color = palette[idx % palette.Length];
            sp.softBodyMaterial = mat;
        }
    }
}
