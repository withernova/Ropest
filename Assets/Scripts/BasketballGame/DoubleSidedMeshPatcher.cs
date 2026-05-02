using System.Collections;
using UnityEngine;

namespace BasketballGame
{
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
