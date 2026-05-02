using System.Collections.Generic;
using UnityEngine;

namespace BasketballGame
{
    public class BasketballSceneBootstrap : MonoBehaviour
    {
        [Header("材质")]
        public Material ballMaterial;
        public Material clothMaterial;
        public Material groundMaterial;

        [Header("球参数")]
        [Range(0.2f, 2f)] public float ballRadius = 0.4f;
        [Range(1, 20)] public int ballResolution = 3;

        [Header("布料参数")]
        [Range(1f, 10f)] public float clothLength = 4f;
        [Range(1f, 10f)] public float clothWidth = 4f;
        [Range(4, 30)] public int clothSegments = 10;
        [Range(4, 30)] public int clothSubdivision = 10;

        [Header("并行求解开关（调试用）")]
        [Tooltip("是否为布料启用图着色(Graph Coloring)并行距离约束求解。\n" +
                 "开启：距离约束按颜色分组并行调度 IJobParallelFor；\n" +
                 "关闭：回退为串行 IJob。")]
        public bool clothUseGraphColoring = false;

        [Tooltip("是否为球体软体启用图着色(Graph Coloring)并行求解（距离约束 + 体积约束）。\n" +
                 "开启：按颜色分组并行；关闭：回退为串行 IJob。")]
        public bool ballUseGraphColoring = false;

        void Awake()
        {
            BuildScene();
        }

        private void BuildScene()
        {
            if (FindObjectOfType<Light>() == null)
            {
                var lightGo = new GameObject("Sun");
                lightGo.transform.rotation = Quaternion.Euler(50f, -30f, 0f);
                var light = lightGo.AddComponent<Light>();
                light.type = LightType.Directional;
                light.intensity = 1.1f;
                light.shadows = LightShadows.Soft;
            }

            // 确保默认材质存在
            var defaultShader = Shader.Find("Universal Render Pipeline/Lit") ?? Shader.Find("Standard");
            if (ballMaterial == null) ballMaterial = new Material(defaultShader) { color = new Color(1f, 0.55f, 0.2f) };
            if (clothMaterial == null) clothMaterial = new Material(defaultShader) { color = new Color(0.6f, 0.85f, 0.4f) };
            if (groundMaterial == null) groundMaterial = new Material(defaultShader) { color = new Color(0.35f, 0.35f, 0.4f) };

            // 布料材质强制双面渲染（URP Lit / Standard 都兼容）
            MakeMaterialDoubleSided(clothMaterial);

            var ground = GameObject.CreatePrimitive(PrimitiveType.Cube);
            ground.name = "Ground";
            ground.transform.position = new Vector3(0, -0.5f, 10f);
            ground.transform.localScale = new Vector3(40f, 1f, 40f);
            ground.GetComponent<MeshRenderer>().sharedMaterial = groundMaterial;
            var groundColliderSrc = ground.AddComponent<AnalyticalColliderSource>();
            groundColliderSrc.colliderType = AnalyticalColliderType.Box;
            groundColliderSrc.autoDetect = true;

            var playerStandGo = new GameObject("PlayerStand");
            playerStandGo.transform.position = new Vector3(0, 0, 0);

            var targetAreaGo = new GameObject("TargetAreaCenter");
            targetAreaGo.transform.position = new Vector3(0, 0, 6f);

            var gmGo = new GameObject("BasketballGameManager");
            var gm = gmGo.AddComponent<BasketballGameManager>();
            gm.playerStand = playerStandGo.transform;
            gm.targetAreaCenter = targetAreaGo.transform;

            Camera cam = Camera.main;
            GameObject camGo;
            if (cam == null)
            {
                camGo = new GameObject("Main Camera");
                camGo.tag = "MainCamera";
                cam = camGo.AddComponent<Camera>();
                camGo.AddComponent<AudioListener>();
            }
            else camGo = cam.gameObject;
            cam.clearFlags = CameraClearFlags.Skybox;
            camGo.transform.position = playerStandGo.transform.position + Vector3.up * 1.6f;
            camGo.transform.rotation = Quaternion.LookRotation(Vector3.forward);

            var camCtrl = camGo.GetComponent<BasketballCameraController>();
            if (camCtrl == null) camCtrl = camGo.AddComponent<BasketballCameraController>();

            var shooter = camGo.GetComponent<BasketballShooter>();
            if (shooter == null) shooter = camGo.AddComponent<BasketballShooter>();
            shooter.shootCamera = cam;
            shooter.ballPrefab = BuildBallPrefab();

            var clothSpawnerGo = new GameObject("ClothSpawner");
            var clothSpawner = clothSpawnerGo.AddComponent<BasketballClothSpawner>();
            clothSpawner.clothPrefab = BuildClothPrefab();

            var hudGo = new GameObject("HUD");
            var hud = hudGo.AddComponent<BasketballHUD>();
            hud.shooter = shooter;
        }

        private GameObject BuildBallPrefab()
        {
            var go = new GameObject("BallPrefab_Template");
            go.SetActive(false);
            var spawner = go.AddComponent<SoftBodyRuntimeSpawner>();
            spawner.shape = SoftBodyRuntimeSpawner.SoftBodyShape.Sphere;
            spawner.sizeX = ballRadius * 2f;
            spawner.sizeY = ballRadius * 2f;
            spawner.sizeZ = ballRadius * 2f;
            spawner.resolution = ballResolution;
            spawner.projectToSphereSurface = true;
            spawner.meshOrigin = new Vector3(-ballRadius, -ballRadius, -ballRadius);
            spawner.softBodyMaterial = ballMaterial;

            spawner.numSubSteps = 10; // 子步更多，跨体碰撞更稳定
            spawner.distanceStiffness = 0f;
            spawner.volumeStiffness = 0f;
            spawner.damping = 0.03f;
            // 关键：球粒子碰撞半径 ≈ 球表面粒子间距的一半，让整个球表面被连续碰撞球覆盖，
            // 避免布料粒子从球面"缝隙"穿入，也让高速球更难 tunneling 穿过布料。
            float ballSurfaceSpacing = (ballRadius * 2f) / Mathf.Max(1, ballResolution);
            spawner.collisionRadius = Mathf.Clamp(ballSurfaceSpacing * 0.6f, 0.08f, 0.4f);
            spawner.friction = 0.35f;
            spawner.renderSubdivisionIterations = 1;
            spawner.gravity = new Vector3(0, -9.8f, 0);
            // 透传图着色并行开关（由 BasketballSceneBootstrap 控制，方便调试对比）
            spawner.useGraphColoring = ballUseGraphColoring;

            go.AddComponent<BasketballBall>();

            // 保留为场景中的不激活 GameObject，供 Shooter Instantiate
            go.transform.SetParent(transform, false);
            return go;
        }

        private GameObject BuildClothPrefab()
        {
            var go = new GameObject("ClothPrefab_Template");
            go.SetActive(false);
            var spawner = go.AddComponent<ClothRuntimeSpawner>();
            spawner.length = clothLength;
            spawner.width = clothWidth;
            spawner.segments = clothSegments;
            spawner.subdivision = clothSubdivision;
            spawner.clothMaterial = clothMaterial;
            spawner.numSubSteps = 10; // 子步更多，跨体碰撞更稳定
            spawner.distanceStiffness = 0f;
            spawner.damping = 0.03f;

            // 关键：碰撞半径必须大于等于"相邻粒子间距的一半"，否则球会从两粒子之间的空洞穿过。
            // 取粒子间距（X/Z 方向更大者）的 ~0.6 倍，让相邻碰撞球足够交叠，形成连续"布面"屏障。
            float spacingX = clothLength / Mathf.Max(1, clothSegments);
            float spacingZ = clothWidth / Mathf.Max(1, clothSubdivision);
            float maxSpacing = Mathf.Max(spacingX, spacingZ);
            spawner.collisionRadius = Mathf.Clamp(maxSpacing * 0.6f, 0.08f, 0.5f);

            spawner.friction = 0.4f;
            spawner.renderSubdivisionIterations = 1;
            spawner.gravity = new Vector3(0, -2f, 0); // 小重力，让布料飘一点但不会垮塌
            // 透传图着色并行开关（由 BasketballSceneBootstrap 控制，方便调试对比）
            spawner.useGraphColoring = clothUseGraphColoring;

            // 固定四个角（让布在空中保持平展，像一块蹦床）
            // 网格索引：segments+1 行 × subdivision+1 列，按 i*(subdivision+1)+j 索引
            int cols = spawner.subdivision + 1;
            int rows = spawner.segments + 1;
            spawner.fixedVertices = new List<int>
            {
                0,                            // (0, 0)
                cols - 1,
                (rows - 1) * cols,            // (segments, 0)
                (rows - 1) * cols + (cols - 1)
            };

            go.AddComponent<BasketballClothTarget>();
            go.AddComponent<DoubleSidedMeshPatcher>();
            go.transform.SetParent(transform, false);
            return go;
        }

        private static void MakeMaterialDoubleSided(Material mat)
        {
            if (mat == null) return;
            // URP Lit / Simple Lit / Unlit: 使用 _Cull 属性（0=Off, 1=Front, 2=Back）
            if (mat.HasProperty("_Cull"))
            {
                mat.SetFloat("_Cull", 0f); // Off = 双面
            }
            // URP 在某些版本用 _CullMode
            if (mat.HasProperty("_CullMode"))
            {
                mat.SetFloat("_CullMode", 0f);
            }
        }
    }
}
