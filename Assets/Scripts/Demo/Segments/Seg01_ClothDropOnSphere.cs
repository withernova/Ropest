using UnityEngine;

namespace Ropest.Demo.Segments
{
    /// <summary>
    /// 段 1：布料自由下落，盖住一个静止的球。
    /// 看点：距离约束 + 球碰撞 + 自然形变。
    /// </summary>
    public class Seg01_ClothDropOnSphere : DemoSegmentBase
    {
        public override string Title => "布料下落覆盖球体";
        public override string Subtitle => "距离约束 + 球面碰撞 · XPBD，开启自碰撞，正反不穿模";
        public override float Duration => 13f;
        public override DemoCameraPreset Camera =>
            DemoCameraPreset.OrbitAndZoom(
                lookAt: new Vector3(0f, 1.0f, 0f),
                distA: 5.5f, distB: 5.0f,
                yawA: 25f, yawB: -10f,
                pitch: 12f, fov: 42f);

        protected override void OnEnter()
        {
            // 静止的球
            var ball = GameObject.CreatePrimitive(PrimitiveType.Sphere);
            ball.name = "Seg01_Ball";
            ball.transform.SetParent(Root, false);
            ball.transform.position = new Vector3(0f, 1.0f, 0f);
            ball.transform.localScale = Vector3.one * 1.6f;
            ball.GetComponent<MeshRenderer>().sharedMaterial = Director.BallMat;
            var src = ball.AddComponent<AnalyticalColliderSource>();
            src.colliderType = AnalyticalColliderType.Sphere;
            src.autoDetect = true;

            // 布料：4×4，segments=40，subdivision=40 → 41×41 ≈ 1681 粒子，盖球绰绰有余
            var clothGo = new GameObject("Seg01_Cloth");
            clothGo.transform.SetParent(Root, false);
            var sp = clothGo.AddComponent<ClothRuntimeSpawner>();
            sp.length = 4f;
            sp.width = 4f;
            sp.segments = 40;
            sp.subdivision = 40;
            sp.meshOrigin = new Vector3(-2f, 3.0f, -2f);
            sp.numSubSteps = 8;
            sp.distanceStiffness = 0f;
            sp.damping = 0.03f;
            sp.collisionRadius = 0.15f;
            sp.friction = 0.5f;
            sp.useGraphColoring = true;
            sp.enableSelfCollision = true;
            sp.renderSubdivisionIterations = 1;
            sp.clothMaterial = Director.ClothMat;
            clothGo.AddComponent<ClothRendererEnabler>();
        }
    }
}
