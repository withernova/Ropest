using System.Collections;
using UnityEngine;

namespace Ropest.Demo
{
    /// <summary>
    /// 录屏 Demo 专用：在布料 spawner 完成 Start() 之后，重新把 MeshRenderer 打开。
    ///
    /// 背景：ClothRuntimeSpawner.Start() 末尾会把 MeshRenderer.enabled 设为 false，
    ///       原本依赖 BasketballClothTarget 之类的脚本在粒子布局完成后再开启。
    ///       Demo 场景没有这种脚本，因此布料永远不可见。
    ///
    /// 做法：本组件在自己的 Start() 里启动协程，等若干帧后强制开启 MeshRenderer。
    ///       延迟一两帧足以让 ClothSimulationSystem 跑过第一次 OnUpdate（让顶点已经被同步到 mesh），
    ///       从而避免"水平初始网格闪现"这个 spawner 想规避的问题。
    /// </summary>
    [DisallowMultipleComponent]
    public class ClothRendererEnabler : MonoBehaviour
    {
        [Tooltip("延迟多少帧后启用 MeshRenderer")]
        public int FrameDelay = 2;

        IEnumerator Start()
        {
            for (int i = 0; i < FrameDelay; i++) yield return null;
            var mr = GetComponent<MeshRenderer>();
            if (mr != null) mr.enabled = true;
        }
    }
}
