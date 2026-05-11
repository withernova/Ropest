using UnityEngine;

namespace Ropest.Demo.Segments
{
    /// <summary>
    /// 段 6：软体 vs 软体碰撞（A 方案：静止大球 + 下落小球）。
    /// 看点：两个软体球发生明显挤压变形 → 分离，不穿透。
    ///
    /// 关键参数说明（为什么这里的 CollisionRadius 和 resolution 要比其他段大）：
    ///   - CrossBodyCollisionSystem 用"粒子-粒子"距离 < (r_i + r_j) 判定接触；
    ///   - 如果表面粒子间距 > 2 * CollisionRadius，粒子会从对方粒子间"穿过去"；
    ///   - 所以必须满足：粒子间距 ≈ size/resolution ≤ 2 * CollisionRadius。
    ///   - 本段 size=2.0 / resolution=12 → 粒子间距 ≈ 0.17，CollisionRadius=0.1 刚好覆盖。
    /// </summary>
    public class Seg06_CrossBodyClothOnSoft : DemoSegmentBase
    {
        public override string Title => "软体 vs 软体 相互碰撞";
        public override string Subtitle => "两软体表面粒子接触 → 挤压形变 → 分离（跨体碰撞系统）";
        public override float Duration => 6f;
        public override DemoCameraPreset Camera =>
            DemoCameraPreset.OrbitAndZoom(
                lookAt: new Vector3(0f, 1.2f, 0f),
                distA: 5.0f, distB: 4.5f,
                yawA: 15f, yawB: 60f,
                pitch: 8f, fov: 42f);

        protected override void OnEnter()
        {
            // 大软体球：位置低，近乎静止（等它自己落到地上稳住）
            var bigGo = new GameObject("Seg06_BigSoftBall");
            bigGo.transform.SetParent(Root, false);
            var big = bigGo.AddComponent<SoftBodyRuntimeSpawner>();
            big.shape = SoftBodyRuntimeSpawner.SoftBodyShape.Sphere;
            big.sizeX = big.sizeY = big.sizeZ = 2.0f;
            big.resolution = 12;               // 表面粒子更密
            big.projectToSphereSurface = true;
            big.meshOrigin = new Vector3(-1.0f, 1.5f, -1.0f); // 半径 1，贴近地面
            big.numSubSteps = 12;
            big.distanceStiffness = 0f;
            big.volumeStiffness = 0.0f;
            big.damping = 0.02f;
            big.collisionRadius = 0.08f;        // 必须 ≥ 粒子间距/2
            big.friction = 0.3f;
            big.useGraphColoring = true;
            big.renderSubdivisionIterations = 2;

            var bigMat = new Material(Director.SoftBodyMat.shader);
            bigMat.color = new Color(0.30f, 0.65f, 0.95f, 1f);
            big.softBodyMaterial = bigMat;

            // 小软体球：从高空落下，正对大球中心
            var smallGo = new GameObject("Seg06_SmallSoftBall");
            smallGo.transform.SetParent(Root, false);
            var small = smallGo.AddComponent<SoftBodyRuntimeSpawner>();
            small.shape = SoftBodyRuntimeSpawner.SoftBodyShape.Sphere;
            small.sizeX = small.sizeY = small.sizeZ = 1.4f;
            small.resolution = 10;
            small.projectToSphereSurface = true;
            small.meshOrigin = new Vector3(-0.7f, 5.5f, -0.7f);
            small.numSubSteps = 12;
            small.distanceStiffness = 0f;
            small.volumeStiffness = 0f;
            small.damping = 0.02f;
            small.collisionRadius = 0.05f;
            small.friction = 0.3f;
            small.useGraphColoring = true;
            small.renderSubdivisionIterations = 2;

            var smallMat = new Material(Director.SoftBodyMat.shader);
            smallMat.color = new Color(0.95f, 0.55f, 0.30f, 1f);
            small.softBodyMaterial = smallMat;
        }
    }
}
