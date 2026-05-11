using System.Collections.Generic;
using UnityEngine;

namespace Ropest.Demo.Segments
{
    /// <summary>
    /// 段 3：布料一侧固定，叠加伪风力，呈现"飘旗"效果。
    /// 看点：风力 + 高频形变。
    /// 实现：复用 Seg02 的悬挂结构，但改成单边固定（旗帜形态），并挂上 WindForceDemo。
    /// </summary>
    public class Seg03_ClothWind : DemoSegmentBase
    {
        public override string Title => "布料受风飘动（旗帜）";
        public override string Subtitle => "周期性水平脉冲驱动 · 顶端高度敏感";
        public override float Duration => 10f;
        public override DemoCameraPreset Camera =>
            DemoCameraPreset.OrbitAndZoom(
                lookAt: new Vector3(1.0f, 2.5f, 0f),
                distA: 7f, distB: 6f,
                yawA: 30f, yawB: 80f,
                pitch: 4f, fov: 42f);

        protected override void OnEnter()
        {
            var clothGo = new GameObject("Seg03_FlagCloth");
            clothGo.transform.SetParent(Root, false);
            var sp = clothGo.AddComponent<ClothRuntimeSpawner>();
            sp.length = 4f;          // X 方向（旗杆挂在 X=0 那一列）
            sp.width = 2.5f;         // Z 方向
            sp.segments = 32;        // 沿 X
            sp.subdivision = 20;     // 沿 Z
            sp.meshOrigin = new Vector3(0f, 4.0f, -1.25f);
            sp.numSubSteps = 12;     // 加大子步，提升稳定性
            sp.distanceStiffness = 0.0001f;  // 略加一点，避免被风过度拉长
            sp.damping = 0.0f;      // 阻尼别太高，否则布料"僵"住
            sp.collisionRadius = 0.05f;
            sp.friction = 0.3f;
            sp.useGraphColoring = true;
            sp.enableSelfCollision = false;
            sp.renderSubdivisionIterations = 1;
            sp.clothMaterial = Director.ClothMat;
            clothGo.AddComponent<ClothRendererEnabler>();

            // 固定 i=0 那一整列（旗杆挂点）→ 索引 0..subdivision
            int subN = sp.subdivision + 1;
            var fixedList = new List<int>(subN);
            for (int j = 0; j < subN; j++) fixedList.Add(j);
            sp.fixedVertices = fixedList;

            // 旗杆视觉
            var pole = GameObject.CreatePrimitive(PrimitiveType.Cylinder);
            pole.name = "Seg03_FlagPole";
            pole.transform.SetParent(Root, false);
            pole.transform.position = new Vector3(0f, 2.5f, 0f);
            pole.transform.localScale = new Vector3(0.08f, 2.5f, 0.08f);
            var capsuleCol = pole.GetComponent<Collider>();
            if (capsuleCol != null) Object.Destroy(capsuleCol);
            var poleMat = new Material(Director.BallMat.shader);
            poleMat.color = new Color(0.18f, 0.18f, 0.20f);
            pole.GetComponent<MeshRenderer>().sharedMaterial = poleMat;

            // 挂风力驱动（拖曳模型，有视觉效果的参数）
            var windGo = new GameObject("Seg03_Wind");
            windGo.transform.SetParent(Root, false);
            var wind = windGo.AddComponent<WindForceDemo>();
            wind.Direction = new Vector3(1f, 0f, 0f);
            wind.TargetSpeed = 6.0f;
            wind.ResponseRate = 4.0f;
            wind.YawWiggleDeg = 25f;
            wind.WiggleFrequency = 1.5f;
            wind.HeightBias = 0.3f;
            wind.MaxSpeed = 10f;
            wind.ApplyEveryFrame = true;
        }
    }
}
