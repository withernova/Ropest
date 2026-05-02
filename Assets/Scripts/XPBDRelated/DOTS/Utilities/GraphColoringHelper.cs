using System.Collections.Generic;
using Unity.Mathematics;

/// <summary>
/// XPBD 约束图着色辅助工具。
///
/// 图着色（Graph Coloring）并行思想：
///   - 把每个约束视作"图的节点"，两个约束若共享任意粒子则相连（即存在写冲突）。
///   - 对该冲突图进行顶点着色：颜色数 = 最大度数 + 1 的上界（贪心算法）。
///   - 同一颜色组内的约束两两不共享粒子，可以用 IJobParallelFor 无锁并行求解。
///
/// 使用方式：
///   1) 调用 ColorEdges / ColorTets 得到每个约束的颜色；
///   2) 调用 GroupByColor 把约束数组原地重排成连续的颜色块，并输出每块的 (Start, Count)。
/// </summary>
public static class GraphColoringHelper
{
    /// <summary>
    /// 对 XPBD 边约束做贪心图着色，返回每条边的颜色索引与颜色总数。
    /// </summary>
    /// <param name="edges">边数组（IndexA / IndexB）。</param>
    /// <param name="numParticles">粒子总数，用于构建"粒子 -> 已用颜色集合"的查找表。</param>
    /// <param name="outColors">长度与 edges 相同，存放每条边分配到的颜色索引。</param>
    /// <returns>颜色总数。</returns>
    public static int ColorEdges(XPBDEdge[] edges, int numParticles, out int[] outColors)
    {
        int n = edges.Length;
        outColors = new int[n];
        if (n == 0) return 0;

        // particleLastColorSlot[p] 记录粒子 p 最近一次被哪个颜色占用，用于快速探测冲突。
        // 为了支持 O(1) 冲突检测，我们维护 particleColorStamp[p, color] = edgeIndex 的稀疏结构：
        // 用 Dictionary<int, int>（每粒子一个 int HashSet）保存该粒子已涉及的颜色集合。
        var particleColors = new HashSet<int>[numParticles];
        for (int i = 0; i < numParticles; i++) particleColors[i] = new HashSet<int>();

        int maxColor = 0;
        for (int e = 0; e < n; e++)
        {
            int a = edges[e].IndexA;
            int b = edges[e].IndexB;

            // 从 0 开始挑一个两端点都未使用过的颜色
            int color = 0;
            while (particleColors[a].Contains(color) || particleColors[b].Contains(color))
            {
                color++;
            }

            outColors[e] = color;
            particleColors[a].Add(color);
            particleColors[b].Add(color);

            if (color + 1 > maxColor) maxColor = color + 1;
        }

        return maxColor;
    }

    /// <summary>
    /// 对四面体体积约束做贪心图着色。两个四面体若共享任一粒子即视为冲突。
    /// </summary>
    public static int ColorTets(Tetrahedron[] tets, int numParticles, out int[] outColors)
    {
        int n = tets.Length;
        outColors = new int[n];
        if (n == 0) return 0;

        var particleColors = new HashSet<int>[numParticles];
        for (int i = 0; i < numParticles; i++) particleColors[i] = new HashSet<int>();

        int maxColor = 0;
        for (int t = 0; t < n; t++)
        {
            var tet = tets[t];
            int i0 = tet.I0, i1 = tet.I1, i2 = tet.I2, i3 = tet.I3;

            int color = 0;
            while (particleColors[i0].Contains(color) ||
                   particleColors[i1].Contains(color) ||
                   particleColors[i2].Contains(color) ||
                   particleColors[i3].Contains(color))
            {
                color++;
            }

            outColors[t] = color;
            particleColors[i0].Add(color);
            particleColors[i1].Add(color);
            particleColors[i2].Add(color);
            particleColors[i3].Add(color);

            if (color + 1 > maxColor) maxColor = color + 1;
        }

        return maxColor;
    }

    /// <summary>
    /// 按颜色对数组做稳定重排，使同一颜色的元素在输出中连续。
    /// 同时输出每个颜色区间的 (Start, Count)。
    /// </summary>
    public static void GroupByColor<T>(T[] items, int[] colors, int numColors,
        out T[] sortedItems, out (int Start, int Count)[] ranges)
    {
        int n = items.Length;
        sortedItems = new T[n];
        ranges = new (int, int)[numColors];

        // 先统计每个颜色的元素数
        var counts = new int[numColors];
        for (int i = 0; i < n; i++) counts[colors[i]]++;

        // 计算起点
        int acc = 0;
        var writeHead = new int[numColors];
        for (int c = 0; c < numColors; c++)
        {
            ranges[c] = (acc, counts[c]);
            writeHead[c] = acc;
            acc += counts[c];
        }

        // 放入数据
        for (int i = 0; i < n; i++)
        {
            int c = colors[i];
            sortedItems[writeHead[c]++] = items[i];
        }
    }
}
