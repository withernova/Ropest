using System.Collections.Generic;
using Ropest.Demo.Segments;
using UnityEngine;

namespace Ropest.Demo
{
    /// <summary>
    /// Demo 录屏总导演。
    /// 使用方法：
    ///   1. 新建一个空场景 DemoRecord.unity（File → New Scene → Empty）。
    ///   2. 场景里新建一个空 GameObject，命名 "DemoDirector"，挂上本脚本。
    ///   3. 直接 Play，即可按时间表自动跑完所有段；用 Unity Recorder 录这个 Game 视图即可。
    ///
    /// 它在 Awake 阶段会程序化地搭好整个场景：
    ///   - 一台主相机（带 DemoCameraRig）
    ///   - 一盏方向光 + 深灰背景（满足你选择的 "纯色深灰背景 + 一盏方向光"）
    ///   - 一块大地板（带 BoxCollider + AnalyticalColliderSource，Y=0）
    ///   - 一个 AnalyticalColliderManager 单例
    ///   - 一个 DemoHud
    ///   - 一个 SegmentRoot 空节点，用于挂当前段产生的临时物体
    /// </summary>
    [DisallowMultipleComponent]
    public class DemoDirector : MonoBehaviour
    {
        [Header("自动开始")]
        public bool AutoStart = true;

        [Header("背景与灯光")]
        public Color BackgroundColor = new Color(0.12f, 0.13f, 0.15f, 1f);
        public Color AmbientColor = new Color(0.25f, 0.27f, 0.30f, 1f);
        public Vector3 SunEulerAngles = new Vector3(50f, 30f, 0f);
        public float SunIntensity = 1.1f;

        [Header("地板")]
        public Vector2 GroundSize = new Vector2(40f, 40f);
        public Color GroundColor = new Color(0.22f, 0.24f, 0.27f, 1f);

        [Header("布料/软体材质（留空则自动生成）")]
        public Material ClothMaterial;
        public Material SoftBodyMaterial;
        public Material BallMaterial;

        // 运行时状态
        readonly List<DemoSegmentBase> _segments = new List<DemoSegmentBase>();
        int _currentIndex = -1;
        DemoSegmentBase _current;
        bool _running;

        Camera _mainCamera;
        DemoCameraRig _cameraRig;
        Transform _segmentRoot;

        public Transform SegmentRoot => _segmentRoot;
        public Material ClothMat => ClothMaterial;
        public Material SoftBodyMat => SoftBodyMaterial;
        public Material BallMat => BallMaterial;

        public string CurrentSegmentTitle => _current != null ? _current.Title : "";
        public string CurrentSegmentSubtitle
        {
            get
            {
                if (_current == null) return "";
                return _isInVisualizationPhase
                    ? "粒子位置 + 速度向量（隐藏 Mesh）"
                    : _current.Subtitle;
            }
        }
        public int CurrentSegmentIndex => _currentIndex;
        public int TotalSegments => _segments.Count;
        public float CurrentSegmentRemaining => _current != null
            ? Mathf.Max(0f, _current.TotalDuration - _segElapsed) : 0f;
        public float CurrentSegmentProgress => _current != null
            ? Mathf.Clamp01(_segElapsed / Mathf.Max(0.001f, _current.TotalDuration)) : 0f;

        float _segElapsed;
        bool _isInVisualizationPhase;
        ParticleVisualizer _activeVisualizer;

        void Awake()
        {
            BuildEnvironment();
            BuildSegments();
        }

        void Start()
        {
            if (AutoStart) BeginRun();
        }

        void BuildEnvironment()
        {
            // ===== 灯光 / 环境 =====
            RenderSettings.ambientMode = UnityEngine.Rendering.AmbientMode.Flat;
            RenderSettings.ambientLight = AmbientColor;
            RenderSettings.fog = false;

            var sunGo = new GameObject("Demo_Sun");
            sunGo.transform.SetParent(transform, false);
            sunGo.transform.eulerAngles = SunEulerAngles;
            var sun = sunGo.AddComponent<Light>();
            sun.type = LightType.Directional;
            sun.intensity = SunIntensity;
            sun.color = Color.white;
            sun.shadows = LightShadows.Soft;
            sun.shadowStrength = 0.6f;

            // ===== 主相机 =====
            var camGo = new GameObject("Demo_MainCamera");
            camGo.transform.SetParent(transform, false);
            _mainCamera = camGo.AddComponent<Camera>();
            _mainCamera.clearFlags = CameraClearFlags.SolidColor;
            _mainCamera.backgroundColor = BackgroundColor;
            _mainCamera.tag = "MainCamera";
            _mainCamera.nearClipPlane = 0.05f;
            _mainCamera.farClipPlane = 200f;

            camGo.AddComponent<AudioListener>();
            _cameraRig = camGo.AddComponent<DemoCameraRig>();
            _cameraRig.TargetCamera = _mainCamera;

            // ===== 地板（视觉 + 解析碰撞体）=====
            var ground = GameObject.CreatePrimitive(PrimitiveType.Cube);
            ground.name = "Demo_Ground";
            ground.transform.SetParent(transform, false);
            ground.transform.position = new Vector3(0f, -0.5f, 0f);
            ground.transform.localScale = new Vector3(GroundSize.x, 1f, GroundSize.y);
            // 删除 Unity Physics 碰撞响应（我们用 AnalyticalColliderSource 给 XPBD 用），但保留 BoxCollider 给自动检测
            var groundCollider = ground.GetComponent<BoxCollider>();
            // 让 XPBD 把它当 Box 解析碰撞
            var src = ground.AddComponent<AnalyticalColliderSource>();
            src.colliderType = AnalyticalColliderType.Box;
            src.autoDetect = true;

            var groundMat = new Material(Shader.Find("Universal Render Pipeline/Lit") ?? Shader.Find("Standard"));
            groundMat.color = GroundColor;
            ground.GetComponent<MeshRenderer>().sharedMaterial = groundMat;

            // ===== AnalyticalColliderManager 单例 =====
            if (AnalyticalColliderManager.Instance == null)
            {
                var mgrGo = new GameObject("AnalyticalColliderManager");
                mgrGo.transform.SetParent(transform, false);
                mgrGo.AddComponent<AnalyticalColliderManager>();
            }

            // ===== 段根节点 =====
            var rootGo = new GameObject("Demo_SegmentRoot");
            rootGo.transform.SetParent(transform, false);
            _segmentRoot = rootGo.transform;

            // ===== HUD =====
            var hudGo = new GameObject("Demo_HUD");
            hudGo.transform.SetParent(transform, false);
            var hud = hudGo.AddComponent<DemoHud>();
            hud.Director = this;

            // ===== 默认材质（如果用户没指定）=====
            if (ClothMaterial == null)
            {
                ClothMaterial = new Material(groundMat.shader);
                ClothMaterial.color = new Color(0.85f, 0.30f, 0.35f, 1f);
                ClothMaterial.name = "Demo_ClothMat";
                // 双面渲染（如果是 URP Lit，开 _Cull = 0）
                if (ClothMaterial.HasProperty("_Cull")) ClothMaterial.SetFloat("_Cull", 0f);
                if (ClothMaterial.HasProperty("_CullMode")) ClothMaterial.SetFloat("_CullMode", 0f);
            }
            if (SoftBodyMaterial == null)
            {
                SoftBodyMaterial = new Material(groundMat.shader);
                SoftBodyMaterial.color = new Color(0.30f, 0.65f, 0.95f, 1f);
                SoftBodyMaterial.name = "Demo_SoftBodyMat";
            }
            if (BallMaterial == null)
            {
                BallMaterial = new Material(groundMat.shader);
                BallMaterial.color = new Color(0.95f, 0.78f, 0.20f, 1f);
                BallMaterial.name = "Demo_BallMat";
            }
        }

        void BuildSegments()
        {
            // 7 段，每段 10s，总 70s。
            // 段内具体参数都在各自类里写死，方便维护 + 阅读。
            _segments.Add(new Seg01_ClothDropOnSphere());
            _segments.Add(new Seg02_ClothHangSwing());
            _segments.Add(new Seg03_ClothWind());
            _segments.Add(new Seg04_SoftBodyCompress());
            _segments.Add(new Seg05_SoftBodyCollision());
            _segments.Add(new Seg06_CrossBodyClothOnSoft());
            _segments.Add(new Seg07_LargeScaleStress());
        }

        public void BeginRun()
        {
            if (_running) return;
            _running = true;
            _currentIndex = -1;
            NextSegment();
        }

        void NextSegment()
        {
            // 退出当前段
            if (_current != null)
            {
                // 如果还在可视化阶段，先把可视化器关掉，恢复 Mesh 渲染
                StopVisualizationIfActive();
                _current.TickExit();
                ClearSegmentRoot();
                _current = null;
            }

            _currentIndex++;
            if (_currentIndex >= _segments.Count)
            {
                _running = false;
                _currentIndex = _segments.Count; // 防越界
                Debug.Log("[DemoDirector] 全部段已演完。");
                return;
            }

            _current = _segments[_currentIndex];
            _current.Bind(this, _segmentRoot);
            _segElapsed = 0f;
            _isInVisualizationPhase = false;
            // 渲染阶段：相机用原 Camera 预设，插值时长 = Duration
            _cameraRig.ApplyPreset(_current.Camera, _current.Duration);
            _current.TickEnter();
            Debug.Log($"[DemoDirector] 段 {_currentIndex + 1}/{_segments.Count}: {_current.Title}" +
                      $"  (Duration={_current.Duration:F1}s + Vis={_current.VisualizationDuration:F1}s)");
        }

        void Update()
        {
            if (!_running || _current == null) return;
            float dt = Time.deltaTime;
            _segElapsed += dt;
            _current.TickUpdate(dt);

            // 检查是否需要进入可视化阶段（仅当 VisualizationDuration > 0 才进入）
            if (!_isInVisualizationPhase
                && _segElapsed >= _current.Duration
                && _current.VisualizationDuration > 0f)
            {
                StartVisualization();
            }

            // 段总时长结束 → 切下一段
            if (_segElapsed >= _current.TotalDuration)
            {
                NextSegment();
            }
        }

        /// <summary>
        /// 进入可视化阶段：硬切重播。
        ///   1. 销毁渲染阶段产生的所有仿真物体（ClearSegmentRoot）；
        ///   2. 重新挂上 ParticleVisualizer，让仿真重生时新出现的 Mesh 立刻被它隐藏；
        ///   3. 用反向相机预设让镜头从原终点反着走回起点；
        ///   4. 重新调段的 TickEnter，用同样参数重生仿真物体（仿真自然从 t=0 开始重新跑）。
        /// </summary>
        void StartVisualization()
        {
            _isInVisualizationPhase = true;

            // 1. 销毁渲染阶段所有物体（含布料/软体/旗杆等所有 Spawner 子物体）
            ClearSegmentRoot();
            _current.TickExit();   // 给段一次清理回调（虽然默认实现是空的，但保险起见）

            // 2. 重挂可视化器（必须在 TickEnter 之前，这样新 Mesh 一出现就被隐藏；
            //    ParticleVisualizer 内部每帧扫一次 MeshRenderer，所以即使 Spawner 在下一帧才 AddComponent
            //    MeshRenderer，下下帧也会被隐藏。）
            var visGo = new GameObject("ParticleVisualizer");
            visGo.transform.SetParent(_segmentRoot, false);
            _activeVisualizer = visGo.AddComponent<ParticleVisualizer>();
            _activeVisualizer.Activate(_mainCamera);

            // 3. 相机反向播放：起止互换，插值时长 = VisualizationDuration
            _cameraRig.ApplyPreset(_current.Camera.Reversed(), _current.VisualizationDuration);

            // 4. 段重新进入：所有 Spawner 用同样参数重生，仿真自动从初始状态开始
            _current.Bind(this, _segmentRoot);
            _current.TickEnter();

            Debug.Log($"[DemoDirector] 段 {_currentIndex + 1} 进入粒子可视化重播（{_current.VisualizationDuration:F1}s，相机反向）");
        }

        void StopVisualizationIfActive()
        {
            if (_activeVisualizer != null)
            {
                _activeVisualizer.Deactivate();
                _activeVisualizer = null;
            }
            _isInVisualizationPhase = false;
        }

        void ClearSegmentRoot()
        {
            if (_segmentRoot == null) return;
            for (int i = _segmentRoot.childCount - 1; i >= 0; i--)
            {
                Destroy(_segmentRoot.GetChild(i).gameObject);
            }
        }
    }
}
