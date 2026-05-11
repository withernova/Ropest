using System.Collections.Generic;
using UnityEngine;

namespace Ropest.Demo.Segments
{
    /// <summary>
    /// 段 2：布料顶端两角固定，下方自由摆动。
    /// 看点：约束传播 + 自然下垂的悬链曲线。
    /// </summary>
    public class Seg02_ClothHangSwing : DemoSegmentBase
    {
        public override string Title => "布料悬挂自由摆动";
        public override string Subtitle => "顶端固定，约束传播形成悬链";
        public override float Duration => 12f;
        public override DemoCameraPreset Camera =>
            DemoCameraPreset.Orbit(
                lookAt: new Vector3(0f, 2.0f, 0f),
                distance: 6.5f,
                yawA: 30f, yawB: 60f,
                pitch: 8f, fov: 42f);

        protected override void OnEnter()
        {
            var clothGo = new GameObject("Seg02_Cloth");
            clothGo.transform.SetParent(Root, false);
            var sp = clothGo.AddComponent<ClothRuntimeSpawner>();
            sp.length = 4f;
            sp.width = 4f;
            sp.segments = 40;
            sp.subdivision = 40;
            sp.meshOrigin = new Vector3(-2f, 4.0f, -2f);
            sp.numSubSteps = 8;
            sp.distanceStiffness = 0f;
            sp.damping = 0.0005f;
            sp.collisionRadius = 0.1f;
            sp.friction = 0.4f;
            sp.useGraphColoring = true;
            sp.enableSelfCollision = true;
            sp.renderSubdivisionIterations = 1;
            sp.clothMaterial = Director.ClothMat;
            clothGo.AddComponent<ClothRendererEnabler>();

            // 固定顶端的两个角点。布料顶点编号规则：
            //   for i in [0, segments]:
            //     for j in [0, subdivision]:
            //       index = i * (subdivision+1) + j
            // meshOrigin 处的角是 index=0；顶端那一边对应 i=0 j 任意 → 索引 0..subdivision。
            // 我们固定两个对角角：j=0 和 j=subdivision（顶边的两个角）。
            int subN = sp.subdivision + 1;
            sp.fixedVertices = new List<int>
            {
                0,                  // (i=0, j=0)
                subN - 1            // (i=0, j=subdivision)
            };
        }
    }
}
