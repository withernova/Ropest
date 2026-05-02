using Unity.Burst;
using Unity.Collections;
using Unity.Jobs;
using Unity.Mathematics;


// 计算每个 body 的粒子质心（表面粒子优先）以及"有效半径"（表面粒子到质心的最大距离）。
// 单线程 IJob：body 数一般很少（< 20），O(N) 累加。
[BurstCompile]
public struct ComputeBodyCentersJob : IJob
{
    [ReadOnly] public NativeArray<float3> GlobalPositions;
    [ReadOnly] public NativeArray<int> GlobalBodyIds;
    [ReadOnly] public NativeArray<byte> GlobalIsSurface;
    public int TotalParticles;

    public NativeArray<float3> BodyCenters;
    // 每个 body 的"有效外接半径"：表面粒子到质心距离的最大值。
    // 用于 body 级 shell 约束：若 |posI - otherCenter| < otherRadius + rI，则说明 i 被 other 的外壳包住，
    // 需要被推到 (otherCenter + normalize(posI - otherCenter) * (otherRadius + rI))
    public NativeArray<float> BodyRadii;

    public void Execute()
    {
        int numBodies = BodyCenters.Length;

        var sumSurf = new NativeArray<float3>(numBodies, Allocator.Temp);
        var cntSurf = new NativeArray<int>(numBodies, Allocator.Temp);
        var sumAll = new NativeArray<float3>(numBodies, Allocator.Temp);
        var cntAll = new NativeArray<int>(numBodies, Allocator.Temp);

        // 第一遍：累加
        for (int i = 0; i < TotalParticles; i++)
        {
            int b = GlobalBodyIds[i];
            if (b < 0 || b >= numBodies) continue;
            float3 p = GlobalPositions[i];
            sumAll[b] += p;
            cntAll[b]++;
            if (GlobalIsSurface[i] != 0)
            {
                sumSurf[b] += p;
                cntSurf[b]++;
            }
        }

        for (int b = 0; b < numBodies; b++)
        {
            if (cntSurf[b] > 0) BodyCenters[b] = sumSurf[b] / cntSurf[b];
            else if (cntAll[b] > 0) BodyCenters[b] = sumAll[b] / cntAll[b];
            else BodyCenters[b] = float3.zero;

            BodyRadii[b] = 0f;
        }

        // 第二遍：算每个 body 的最大表面粒子到质心距离
        for (int i = 0; i < TotalParticles; i++)
        {
            if (GlobalIsSurface[i] == 0) continue;
            int b = GlobalBodyIds[i];
            if (b < 0 || b >= numBodies) continue;
            float d = math.distance(GlobalPositions[i], BodyCenters[b]);
            if (d > BodyRadii[b]) BodyRadii[b] = d;
        }

        sumSurf.Dispose();
        cntSurf.Dispose();
        sumAll.Dispose();
        cntAll.Dispose();
    }
}

[BurstCompile]
public struct BuildCrossBodyHashJob : IJobParallelFor
{
    [ReadOnly] public NativeArray<float3> GlobalPositions;
    // 所有粒子都写入哈希（包括软体内部粒子），保证穿透总能被检测到
    [ReadOnly] public NativeArray<byte> GlobalIsSurface;
    public NativeParallelMultiHashMap<int, int>.ParallelWriter HashMap;
    public float CellSize;

    public void Execute(int i)
    {
        float invCell = 1f / CellSize;
        int hash = SpatialHash(GlobalPositions[i], invCell);
        HashMap.Add(hash, i);
    }

    public static int SpatialHash(float3 pos, float invCell)
    {
        int x = (int)math.floor(pos.x * invCell);
        int y = (int)math.floor(pos.y * invCell);
        int z = (int)math.floor(pos.z * invCell);
        return x * 73856093 ^ y * 19349663 ^ z * 83492791;
    }
}

[BurstCompile]
public struct CrossBodyCollisionResolveJob : IJobParallelFor
{
    [ReadOnly] public NativeArray<float3> GlobalPositions;
    [ReadOnly] public NativeArray<float3> GlobalPrevPositions;
    [ReadOnly] public NativeArray<float> GlobalInvMasses;
    [ReadOnly] public NativeArray<float> GlobalRadii;
    [ReadOnly] public NativeArray<int> GlobalBodyIds;
    [ReadOnly] public NativeArray<byte> GlobalIsSurface;
    [ReadOnly] public NativeArray<float3> BodyCenters;
    // 每个 body 的有效半径（表面粒子最大到质心距离）
    [ReadOnly] public NativeArray<float> BodyRadii;
    // 每个 body 是否是布料（1=布料，0=软体）
    // 涉及布料的碰撞不能用 body 级 shell 约束（布料是平面/薄片，外接球半径很大且质心与法线无关），
    // 也不能用"posI - opponentCenter"作方向（方向可能完全错，导致粒子被沿布料平面弹飞）。
    // 对这种情况，回退到经典的"posI - posJ"粒子-粒子连线方向（贴近布料局部法线）。
    [ReadOnly] public NativeArray<byte> BodyIsCloth;
    [ReadOnly] public NativeParallelMultiHashMap<int, int> HashMap;

    public NativeArray<float3> GlobalCorrections;

    public float CellSize;

    // 单次迭代单个粒子的最大修正幅度 = MaxCorrectionRatio * (rI + rJ)
    public float MaxCorrectionRatio;

    // === 核心思路：body 级 shell 约束 + 粒子级 overlap，取两者较大修正 ===
    //
    // 问题回顾：之前只用粒子-粒子 overlap 求解，
    // 球 A 半个球穿入 B 时，A 的表面粒子 i 的邻域里全是 B 的内部粒子，
    // 粒子级 overlap 最多 ~minDist=rI+rJ（例如 0.32），但要把 i 推出 B 需要走 ~球半径 0.4 米，
    // 8 次迭代的上限也才 ~8*0.5*0.32 = 1.28 米，虽然理论够，但每次要搜到最深邻居且迭代是渐进的，
    // 实际很难把深穿透解开。
    //
    // 新增：body 级 shell 约束
    //   - 只要 i 是对方 body 的"外接球"内点（|posI - otherCenter| < otherRadius + rI），
    //     就直接把 i 推到外接球表面：newPos = otherCenter + normalize(posI - otherCenter) * (otherRadius + rI)
    //   - 这个约束对"半球穿入"尤其有效——它一次就能把 i 推到 B 真正的外缘，不再依赖邻居搜索。
    //   - 对浅接触（粒子级 overlap 已足够）也无害：粒子级修正可能更大，取较大者。
    //
    // 为什么不会破坏形变：shell 约束只把 i 推到"对方的最外表面"，之后靠本身的距离/体积约束恢复形状。

    public void Execute(int i)
    {
        if (GlobalIsSurface[i] == 0)
        {
            GlobalCorrections[i] = float3.zero;
            return;
        }

        float wI = GlobalInvMasses[i];
        if (wI <= 0f)
        {
            GlobalCorrections[i] = float3.zero;
            return;
        }

        float3 posI = GlobalPositions[i];
        float3 prevPosI = GlobalPrevPositions[i];
        int bodyI = GlobalBodyIds[i];
        float rI = GlobalRadii[i];

        float invCell = 1f / CellSize;
        int cx = (int)math.floor(posI.x * invCell);
        int cy = (int)math.floor(posI.y * invCell);
        int cz = (int)math.floor(posI.z * invCell);

        // === 遍历邻域，找最深"粒子级 overlap"，同时记录所有遇到的对方 body（用于 shell 约束）===
        float maxOverlap = 0f;
        float maxOverlapMinDist = 0f;
        int deepestJ = -1;
        int deepestBody = -1;
        float shareIForMax = 0.5f;

        // 用位图/简单数组记录邻域内遇到的对方 body（body 数最多 ~32 足够）
        // 记录前 4 个对方 body 即可，罕见情况下 i 同时接触 >4 个 body 时只处理前 4 个
        const int kMaxNearbyBodies = 4;
        int nearbyCount = 0;
        var nearbyBodies = new int4(-1, -1, -1, -1);

        for (int dx = -1; dx <= 1; dx++)
        {
            for (int dy = -1; dy <= 1; dy++)
            {
                for (int dz = -1; dz <= 1; dz++)
                {
                    int hash = (cx + dx) * 73856093 ^ (cy + dy) * 19349663 ^ (cz + dz) * 83492791;
                    if (!HashMap.TryGetFirstValue(hash, out int j, out var it)) continue;
                    do
                    {
                        if (j == i) continue;
                        int bodyJ = GlobalBodyIds[j];
                        if (bodyJ == bodyI) continue;

                        float wJ = GlobalInvMasses[j];
                        if (wI + wJ < 1e-8f) continue;

                        float shareI = (wJ > 0f) ? 0.5f : 1.0f;

                        float rJ = GlobalRadii[j];
                        float minDist = rI + rJ;
                        float minDistSq = minDist * minDist;

                        float3 diff = posI - GlobalPositions[j];
                        float distSq = math.lengthsq(diff);
                        if (distSq >= minDistSq) continue;

                        // 记录邻域对方 body
                        bool known =
                            (nearbyBodies.x == bodyJ) ||
                            (nearbyBodies.y == bodyJ) ||
                            (nearbyBodies.z == bodyJ) ||
                            (nearbyBodies.w == bodyJ);
                        if (!known && nearbyCount < kMaxNearbyBodies)
                        {
                            if (nearbyCount == 0) nearbyBodies.x = bodyJ;
                            else if (nearbyCount == 1) nearbyBodies.y = bodyJ;
                            else if (nearbyCount == 2) nearbyBodies.z = bodyJ;
                            else nearbyBodies.w = bodyJ;
                            nearbyCount++;
                        }

                        float dist = (distSq < 1e-12f) ? 0f : math.sqrt(distSq);
                        float overlap = minDist - dist;

                        if (overlap > maxOverlap)
                        {
                            maxOverlap = overlap;
                            maxOverlapMinDist = minDist;
                            shareIForMax = shareI;
                            deepestJ = j;
                            deepestBody = bodyJ;
                        }
                    } while (HashMap.TryGetNextValue(out j, ref it));
                }
            }
        }

        if (nearbyCount == 0)
        {
            GlobalCorrections[i] = float3.zero;
            return;
        }

        // 自己或对方是布料？布料是平面结构，不能应用 body shell 约束与"质心方向"假设
        bool selfIsCloth = BodyIsCloth[bodyI] != 0;
        bool deepestIsCloth = (deepestBody >= 0) && (BodyIsCloth[deepestBody] != 0);
        bool involvesClothForParticle = selfIsCloth || deepestIsCloth;

        // === 从接触到的每个对方 body 取 shell 约束，选"推开量最大"的那个 ===
        // Shell 约束（仅对"明显穿入 body 内部"的场景生效）：
        //   若 |posI - otherCenter| < otherRadius，说明 i 已在对方 body 的外接球内部，
        //   直接把 i 推到外接球表面 (otherCenter + dir * otherRadius)。
        //   推开量 shellPush = otherRadius - |posI - otherCenter|。
        //
        // 关键：若自己是布料，或候选对方 body 是布料，shell 完全跳过——
        //   布料的"外接球"是沿平面铺开的大球，对其做 shell push 会把粒子沿平面推飞到远处。
        float bestShellPush = 0f;
        float3 bestShellDir = float3.zero;
        int bestShellBody = -1;

        if (!selfIsCloth)
        {
            for (int k = 0; k < nearbyCount; k++)
            {
                int b;
                if (k == 0) b = nearbyBodies.x;
                else if (k == 1) b = nearbyBodies.y;
                else if (k == 2) b = nearbyBodies.z;
                else b = nearbyBodies.w;
                if (b < 0) continue;
                if (BodyIsCloth[b] != 0) continue; // 对方是布料 → 不做 shell 推开

                float3 oCenter = BodyCenters[b];
                float oRadius = BodyRadii[b];
                float targetDist = oRadius; // i 应该在对方 body 外接球之外

                float3 dCenter = posI - oCenter;
                float dLenSq = math.lengthsq(dCenter);
                if (dLenSq >= targetDist * targetDist) continue; // 未穿入 shell

                float dLen = (dLenSq < 1e-12f) ? 0f : math.sqrt(dLenSq);
                float shellPush = targetDist - dLen;
                if (shellPush > bestShellPush)
                {
                    bestShellPush = shellPush;
                    bestShellBody = b;
                    if (dLen > 1e-6f)
                    {
                        bestShellDir = dCenter / dLen;
                    }
                    else
                    {
                        // i 恰好在对方质心上：用前一帧方向兜底，否则伪随机
                        bestShellDir = float3.zero;
                    }
                }
            }
        }

        // === 选择最终的修正量：body shell 修正通常比粒子级大（半球穿透时尤其明显），
        //     粒子级修正在浅接触时更精确，取两者较大者沿对应方向推开 ===
        float3 finalCorr = float3.zero;

        // 粒子级修正：
        //   - 双方都是软体：方向用 "对方质心 → i"，鲁棒的"远离对方"方向
        //   - 任一方是布料：方向用 "posJ_deepest → i"，即粒子-粒子连线，贴近布料局部法线
        if (maxOverlap > 0f && deepestBody >= 0)
        {
            float3 dir;
            if (involvesClothForParticle)
            {
                // 布料碰撞：粒子-粒子连线方向（更接近布料局部法线）
                float3 dDiff = posI - GlobalPositions[deepestJ];
                float dDiffLenSq = math.lengthsq(dDiff);
                if (dDiffLenSq > 1e-10f)
                {
                    dir = dDiff * math.rsqrt(dDiffLenSq);
                }
                else
                {
                    // 距离为零的退化：用前一帧位置差做参考
                    float3 refDiff = prevPosI - GlobalPrevPositions[deepestJ];
                    float refLenSq = math.lengthsq(refDiff);
                    if (refLenSq > 1e-12f) dir = refDiff * math.rsqrt(refLenSq);
                    else
                    {
                        int seed = (i - deepestJ);
                        dir = math.normalize(new float3(
                            ((seed * 12.9898f) % 1f) - 0.5f,
                            ((seed * 78.2330f) % 1f) - 0.5f + 0.1f,
                            ((seed * 37.7190f) % 1f) - 0.5f));
                    }
                }
            }
            else
            {
                // 双方软体：用对方质心方向（稳健）
                float3 oCenter = BodyCenters[deepestBody];
                float3 dCenter = posI - oCenter;
                float dLenSq = math.lengthsq(dCenter);
                if (dLenSq > 1e-10f)
                {
                    dir = dCenter * math.rsqrt(dLenSq);
                }
                else
                {
                    float3 refDiff = prevPosI - GlobalPrevPositions[deepestJ];
                    float refLenSq = math.lengthsq(refDiff);
                    if (refLenSq > 1e-12f) dir = refDiff * math.rsqrt(refLenSq);
                    else
                    {
                        int seed = (i - deepestJ);
                        dir = math.normalize(new float3(
                            ((seed * 12.9898f) % 1f) - 0.5f,
                            ((seed * 78.2330f) % 1f) - 0.5f + 0.1f,
                            ((seed * 37.7190f) % 1f) - 0.5f));
                    }
                }
            }
            float partMag = shareIForMax * maxOverlap;
            float partCap = MaxCorrectionRatio * maxOverlapMinDist;
            if (partMag > partCap) partMag = partCap;
            finalCorr = dir * partMag;
        }

        // body shell 修正：仅在不涉及布料时启用（前面已在收集 shell 候选时过滤）
        if (bestShellBody >= 0 && bestShellPush > math.length(finalCorr))
        {
            float3 dir = bestShellDir;
            if (math.lengthsq(dir) < 1e-10f)
            {
                // 方向退化：用 prev 位置兜底
                float3 refDir = prevPosI - BodyCenters[bestShellBody];
                float rl = math.lengthsq(refDir);
                if (rl > 1e-12f) dir = refDir * math.rsqrt(rl);
                else dir = new float3(0, 1, 0);
            }
            float shellMag = 0.5f * bestShellPush;
            float shellCap = 2.0f * (rI * 2f);
            if (shellMag > shellCap) shellMag = shellCap;
            finalCorr = dir * shellMag;
        }

        GlobalCorrections[i] = finalCorr;
    }
}

[BurstCompile]
public struct ApplyGlobalCorrectionJob : IJobParallelFor
{
    public NativeArray<float3> GlobalPositions;
    public NativeArray<float3> GlobalCorrections;

    public void Execute(int i)
    {
        GlobalPositions[i] += GlobalCorrections[i];
        GlobalCorrections[i] = float3.zero;
    }
}
