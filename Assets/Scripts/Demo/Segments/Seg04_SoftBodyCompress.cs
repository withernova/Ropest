using UnityEngine;

namespace Ropest.Demo.Segments
{
    /// <summary>
    /// 段 4：软体球从空中下落，砸到地面后压缩、回弹。
    /// 看点：体积约束 + 距离约束 + 地面碰撞 → 弹性形变。
    /// </summary>
    public class Seg04_SoftBodyCompress : DemoSegmentBase
    {
        public override string Title => "软体下落压缩与回弹";
        public override string Subtitle => "距离约束 + 体积约束（XPBD 体积守恒）";
        public override float Duration => 4f;
        public override DemoCameraPreset Camera =>
            DemoCameraPreset.OrbitAndZoom(
                lookAt: new Vector3(0f, 0.8f, 0f),
                distA: 5.0f, distB: 4.5f,
                yawA: -25f, yawB: 20f,
                pitch: 8f, fov: 42f);

        protected override void OnEnter()
        {
            var go = new GameObject("Seg04_SoftBall");
            go.transform.SetParent(Root, false);
            var sp = go.AddComponent<SoftBodyRuntimeSpawner>();
            sp.shape = SoftBodyRuntimeSpawner.SoftBodyShape.Sphere;
            sp.sizeX = 1.6f;
            sp.sizeY = 1.6f;
            sp.sizeZ = 1.6f;
            sp.resolution = 8;
            sp.projectToSphereSurface = true;
            sp.meshOrigin = new Vector3(-0.8f, 4.5f, -0.8f);
            sp.numSubSteps = 12;
            sp.distanceStiffness = 0.0005f;
            sp.volumeStiffness = 0.0005f;
            sp.damping = 0.02f;
            sp.collisionRadius = 0.05f;
            sp.friction = 0.4f;
            sp.useGraphColoring = true;
            sp.renderSubdivisionIterations = 1;
            sp.softBodyMaterial = Director.SoftBodyMat;
        }
    }
}
