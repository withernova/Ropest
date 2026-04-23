using System.Collections;
using UnityEngine;

namespace BasketballGame
{
    /// <summary>
    /// 把 MeshFilter 的 Mesh 改造为"双面可见"：
    /// 做法是追加一份反向三角形（交换两个索引顺序）到同一 SubMesh 中。
    /// 
    /// 使用方式：挂在布料 Prefab 上即可，本脚本会等 ClothRuntimeSpawner 把 mesh 建好后执行一次，
    /// 与具体材质/shader 无关，比设置 _Cull 更稳妥。
    ///
    /// 注意：
    /// 1. 这只改变"渲染用"的三角形索引，不影响 XPBD 仿真（仿真在自己的 simMesh 上进行）
    /// 2. MeshUpdateSystem 会每帧 SetVertices，但不会重新 SetTriangles，所以我们 patch 一次即可
    /// 3. 不要加 [RequireComponent(typeof(MeshFilter))]，因为 ClothRuntimeSpawner 会在 Start 中
    ///    动态 AddComponent<MeshFilter>，若 GameObject 已有 MeshFilter 会返回 null 导致 NRE。
    /// </summary>
    public class DoubleSidedMeshPatcher : MonoBehaviour
    {
        [Tooltip("等待几帧后再 patch（给 ClothRuntimeSpawner 的 Start 留时间）")]
        public int waitFrames = 1;

        private bool _patched;

        void Start()
        {
            StartCoroutine(PatchAfterFrames());
        }

        private IEnumerator PatchAfterFrames()
        {
            for (int i = 0; i < waitFrames; i++) yield return null;
            PatchNow();
        }

        public void PatchNow()
        {
            if (_patched) return;
            var mf = GetComponent<MeshFilter>();
            if (mf == null) return;
            // 必须用 sharedMesh：用 .mesh 会触发 Unity 运行时克隆一份新 mesh 实例给 MeshFilter，
            // 但 ClothMeshUpdateSystem 用的是 ManagedMeshReference.Mesh（原始 clothMesh 引用）
            // 若在此处克隆，MeshFilter 显示的是克隆体（我们 patch 的那份），
            // 而 XPBD 每帧写入的是原始 clothMesh，导致布料"看起来从不更新顶点"
            var mesh = mf.sharedMesh;
            if (mesh == null) return;

            int subMeshCount = mesh.subMeshCount;
            if (subMeshCount == 0) return;

            // 将每个 SubMesh 的三角形追加一份反向
            for (int s = 0; s < subMeshCount; s++)
            {
                var tris = mesh.GetTriangles(s);
                if (tris == null || tris.Length == 0) continue;
                int origLen = tris.Length;
                var newTris = new int[origLen * 2];
                System.Array.Copy(tris, 0, newTris, 0, origLen);
                for (int i = 0; i < origLen; i += 3)
                {
                    // 反向绕序：i0, i2, i1
                    newTris[origLen + i] = tris[i];
                    newTris[origLen + i + 1] = tris[i + 2];
                    newTris[origLen + i + 2] = tris[i + 1];
                }
                mesh.SetTriangles(newTris, s, false);
            }

            // 重新计算法线以避免背面光照异常（加法同位置 + 反向法线会抵消，
            // 这里保留原法线即可，背面会使用同一法线导致"暗面"，这是可接受的效果）
            _patched = true;
        }
    }
}
