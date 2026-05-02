using System.Collections.Generic;
using UnityEngine;

public static class SubdivisionMeshHelper
{
    public struct BindingData
    {
        public int SimI0, SimI1, SimI2; // 原始模拟三角形的三个顶点索引
        public float U, V, W;           // 相对于原始三角形的重心坐标
    }

    public static void SubdivideWithBindings(
        Vector3[] srcVertices, int[] srcTriangles, Vector2[] srcUVs, int iterations,
        out Vector3[] outVertices, out int[] outTriangles, out Vector2[] outUVs,
        out BindingData[] outBindings)
    {
        // === 第一步：纯几何细分（只管顶点位置、UV、三角形索引） ===
        var verts = new List<Vector3>(srcVertices);
        var tris = new List<int>(srcTriangles);
        var uvs = new List<Vector2>();
        bool hasUV = srcUVs != null && srcUVs.Length == srcVertices.Length;
        if (hasUV) uvs.AddRange(srcUVs);

        int numOrigTris = srcTriangles.Length / 3;

        // 每个细分三角形对应的原始三角形索引
        var triOrigIndex = new List<int>(numOrigTris);
        for (int t = 0; t < numOrigTris; t++)
            triOrigIndex.Add(t);

        var triBarys = new List<(Vector3 b0, Vector3 b1, Vector3 b2)>(numOrigTris);
        for (int t = 0; t < numOrigTris; t++)
        {
            // 原始三角形的三个顶点的重心坐标就是 (1,0,0), (0,1,0), (0,0,1)
            triBarys.Add((
                new Vector3(1, 0, 0),
                new Vector3(0, 1, 0),
                new Vector3(0, 0, 1)
            ));
        }

        for (int iter = 0; iter < iterations; iter++)
        {
            var newTris = new List<int>();
            var newTriOrigIndex = new List<int>();
            var newTriBarys = new List<(Vector3, Vector3, Vector3)>();
            // 边中点缓存：(min, max) -> 新顶点索引（只用于几何位置和UV）
            var edgeMidpoints = new Dictionary<(int, int), int>();

            int triCount = tris.Count / 3;
            for (int t = 0; t < triCount; t++)
            {
                int i0 = tris[t * 3 + 0];
                int i1 = tris[t * 3 + 1];
                int i2 = tris[t * 3 + 2];
                int origTri = triOrigIndex[t];
                var (b0, b1, b2) = triBarys[t];

                // 获取或创建三条边的中点（只管几何位置）
                int m01 = GetOrCreateMidpoint(edgeMidpoints, verts, uvs, hasUV, i0, i1);
                int m12 = GetOrCreateMidpoint(edgeMidpoints, verts, uvs, hasUV, i1, i2);
                int m20 = GetOrCreateMidpoint(edgeMidpoints, verts, uvs, hasUV, i2, i0);

                // 中点的重心坐标 = 两端点重心坐标的平均（在同一个原始三角形内，这是精确的）
                Vector3 bm01 = (b0 + b1) * 0.5f;
                Vector3 bm12 = (b1 + b2) * 0.5f;
                Vector3 bm20 = (b2 + b0) * 0.5f;

                // 4个子三角形，都继承同一个原始三角形索引
                // 子三角形0: (i0, m01, m20)
                newTris.AddRange(new[] { i0, m01, m20 });
                newTriOrigIndex.Add(origTri);
                newTriBarys.Add((b0, bm01, bm20));

                // 子三角形1: (m01, i1, m12)
                newTris.AddRange(new[] { m01, i1, m12 });
                newTriOrigIndex.Add(origTri);
                newTriBarys.Add((bm01, b1, bm12));

                // 子三角形2: (m20, m12, i2)
                newTris.AddRange(new[] { m20, m12, i2 });
                newTriOrigIndex.Add(origTri);
                newTriBarys.Add((bm20, bm12, b2));

                // 子三角形3: (m01, m12, m20) — 中心三角形
                newTris.AddRange(new[] { m01, m12, m20 });
                newTriOrigIndex.Add(origTri);
                newTriBarys.Add((bm01, bm12, bm20));
            }

            tris = newTris;
            triOrigIndex = newTriOrigIndex;
            triBarys = newTriBarys;
        }

        outVertices = verts.ToArray();
        outTriangles = tris.ToArray();
        outUVs = hasUV ? uvs.ToArray() : null;

        outBindings = new BindingData[verts.Count];
        var assigned = new bool[verts.Count];

        int finalTriCount = tris.Count / 3;
        for (int t = 0; t < finalTriCount; t++)
        {
            int origTri = triOrigIndex[t];
            int si0 = srcTriangles[origTri * 3 + 0];
            int si1 = srcTriangles[origTri * 3 + 1];
            int si2 = srcTriangles[origTri * 3 + 2];

            var (b0, b1, b2) = triBarys[t];

            int vi0 = tris[t * 3 + 0];
            int vi1 = tris[t * 3 + 1];
            int vi2 = tris[t * 3 + 2];

            if (!assigned[vi0])
            {
                assigned[vi0] = true;
                outBindings[vi0] = new BindingData
                {
                    SimI0 = si0, SimI1 = si1, SimI2 = si2,
                    U = b0.x, V = b0.y, W = b0.z
                };
            }
            if (!assigned[vi1])
            {
                assigned[vi1] = true;
                outBindings[vi1] = new BindingData
                {
                    SimI0 = si0, SimI1 = si1, SimI2 = si2,
                    U = b1.x, V = b1.y, W = b1.z
                };
            }
            if (!assigned[vi2])
            {
                assigned[vi2] = true;
                outBindings[vi2] = new BindingData
                {
                    SimI0 = si0, SimI1 = si1, SimI2 = si2,
                    U = b2.x, V = b2.y, W = b2.z
                };
            }
        }

        // 安全检查：所有顶点都应该被分配了绑定
        for (int i = 0; i < verts.Count; i++)
        {
            if (!assigned[i])
            {
                Debug.LogWarning($"[SubdivisionMeshHelper] 顶点 {i} 未被分配绑定，使用默认值");
                outBindings[i] = new BindingData
                {
                    SimI0 = srcTriangles[0],
                    SimI1 = srcTriangles[1],
                    SimI2 = srcTriangles[2],
                    U = 1, V = 0, W = 0
                };
            }
        }
    }

    static int GetOrCreateMidpoint(
        Dictionary<(int, int), int> cache,
        List<Vector3> verts, List<Vector2> uvs, bool hasUV,
        int a, int b)
    {
        var key = (Mathf.Min(a, b), Mathf.Max(a, b));
        if (cache.TryGetValue(key, out int midIdx))
            return midIdx;

        midIdx = verts.Count;
        verts.Add((verts[a] + verts[b]) * 0.5f);
        if (hasUV)
            uvs.Add((uvs[a] + uvs[b]) * 0.5f);

        cache[key] = midIdx;
        return midIdx;
    }
}
