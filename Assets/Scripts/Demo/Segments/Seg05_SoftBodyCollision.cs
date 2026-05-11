using Unity.Mathematics;
using UnityEngine;

namespace Ropest.Demo.Segments
{
    /// <summary>
    /// 段 5：多个软体（球 + 立方体）从不同高度落下，互相碰撞、碰到地面、产生不同形变。
    /// 看点：几何形状的影响（球 vs 立方体）、不同 stiffness 的对比、跨物体碰撞的稳定性。
    /// </summary>
    public class Seg05_SoftBodyCollision : DemoSegmentBase
    {
        public override string Title => "软体多体碰撞与形变";
        public override string Subtitle => "球与立方体不同几何 · 不同柔度对比";
        public override float Duration => 5f;
        public override DemoCameraPreset Camera =>
            DemoCameraPreset.OrbitAndZoom(
                lookAt: new Vector3(0f, 0.8f, 0f),
                distA: 6.5f, distB: 5.5f,
                yawA: 35f, yawB: -25f,
                pitch: 10f, fov: 45f);

        protected override void OnEnter()
        {
            // 球 1：偏软
            SpawnSoftSphere("Seg05_SoftBall_A",
                position: new Vector3(-1.2f, 4.0f, 0f),
                size: 1.4f,
                resolution: 10,
                stiffness: 0.000001f,
                mat: Director.SoftBodyMat);

            // 球 2：稍硬
            SpawnSoftSphere("Seg05_SoftBall_B",
                position: new Vector3(1.2f, 5.5f, 0.3f),
                size: 1.2f,
                resolution: 10,
                stiffness: 0.0000001f,
                mat: Director.SoftBodyMat);

            // 立方体：橙色
            SpawnSoftCube("Seg05_SoftCube",
                position: new Vector3(0f, 7.5f, -0.3f),
                size: 1.3f,
                resolution: 8,
                stiffness: 0.000005f);
        }

        void SpawnSoftSphere(string name, Vector3 position, float size, int resolution, float stiffness, Material mat)
        {
            var go = new GameObject(name);
            go.transform.SetParent(Root, false);
            var sp = go.AddComponent<SoftBodyRuntimeSpawner>();
            sp.shape = SoftBodyRuntimeSpawner.SoftBodyShape.Sphere;
            sp.sizeX = sp.sizeY = sp.sizeZ = size;
            sp.resolution = resolution;
            sp.projectToSphereSurface = true;
            sp.meshOrigin = position - new Vector3(size, size, size) * 0.5f;
            sp.numSubSteps = 12;
            sp.distanceStiffness = stiffness;
            sp.volumeStiffness = stiffness;
            sp.damping = 0.02f;
            // 必须 ≥ 粒子间距/2（粒子间距 ≈ size/resolution）。否则两球会互相穿透。
            sp.collisionRadius = math.max(0.08f, (size / resolution) * 0.6f);
            sp.friction = 0.4f;
            sp.useGraphColoring = true;
            sp.renderSubdivisionIterations = 1;
            sp.softBodyMaterial = mat;
        }

        void SpawnSoftCube(string name, Vector3 position, float size, int resolution, float stiffness)
        {
            var go = new GameObject(name);
            go.transform.SetParent(Root, false);
            var sp = go.AddComponent<SoftBodyRuntimeSpawner>();
            sp.shape = SoftBodyRuntimeSpawner.SoftBodyShape.Cuboid;
            sp.sizeX = sp.sizeY = sp.sizeZ = size;
            sp.resolution = resolution;
            sp.meshOrigin = position - new Vector3(size, size, size) * 0.5f;
            sp.numSubSteps = 12;
            sp.distanceStiffness = stiffness;
            sp.volumeStiffness = stiffness;
            sp.damping = 0.02f;
            sp.collisionRadius = math.max(0.08f, (size / resolution) * 0.6f);
            sp.friction = 0.4f;
            sp.useGraphColoring = true;
            sp.renderSubdivisionIterations = 1;
            // 立方体用稍微不同的颜色，让画面里能区分
            var cubeMat = new Material(Director.SoftBodyMat.shader);
            cubeMat.color = new Color(0.95f, 0.55f, 0.30f, 1f);
            sp.softBodyMaterial = cubeMat;
        }
    }
}
