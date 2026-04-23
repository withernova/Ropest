using Unity.Burst;
using Unity.Collections;
using Unity.Jobs;
using Unity.Mathematics;
using System.Runtime.InteropServices;

// ============================================================
// 渲染网格插值 Jobs
// 位置：纯线性重心坐标插值（物理正确）
// 法线：用模拟网格的顶点法线进行重心坐标插值（视觉平滑）
// ============================================================

/// <summary>
/// 渲染顶点绑定数据（Burst友好的blittable结构）
/// 内存布局必须与 RenderVertexBinding 完全一致（int,int,int,float,float,float 顺序），
/// 以便直接通过 DynamicBuffer.Reinterpret 零拷贝复用。
/// </summary>
[StructLayout(LayoutKind.Sequential)]
public struct RenderBindingData
{
    public int SimI0;   // 模拟三角形顶点索引0
    public int SimI1;   // 模拟三角形顶点索引1
    public int SimI2;   // 模拟三角形顶点索引2
    public float U;     // 重心坐标u（对应SimI0）
    public float V;     // 重心坐标v（对应SimI1）
    public float W;     // 重心坐标w（对应SimI2）
}

/// <summary>
/// 位置 + 法线 插值Job
/// 
/// 位置：纯线性重心坐标插值
///   renderPos = u*p0 + v*p1 + w*p2
/// 
/// 法线：用模拟网格的顶点法线（面积加权平均，已平滑）进行重心坐标插值
///   renderNormal = normalize(u*n0 + v*n1 + w*n2)
/// 
/// 这样即使位置在平面上，法线也是平滑过渡的，
/// 渲染时光照计算会产生平滑的视觉效果（类似Phong Shading）
/// </summary>
[BurstCompile]
public struct InterpolatePositionAndNormalJob : IJobParallelFor
{
    [ReadOnly] public NativeArray<float3> SimPositions;     // 模拟粒子位置
    [ReadOnly] public NativeArray<float3> SimNormals;       // 模拟粒子法线（面积加权平均）
    [ReadOnly] public NativeArray<RenderBindingData> Bindings;
    [WriteOnly] public NativeArray<float3> RenderPositions; // 输出：渲染顶点位置
    [WriteOnly] public NativeArray<float3> RenderNormals;   // 输出：渲染顶点法线

    public void Execute(int i)
    {
        var b = Bindings[i];
        float3 p0 = SimPositions[b.SimI0];
        float3 p1 = SimPositions[b.SimI1];
        float3 p2 = SimPositions[b.SimI2];

        // 位置：纯线性插值
        RenderPositions[i] = b.U * p0 + b.V * p1 + b.W * p2;

        // 法线：用模拟顶点的平滑法线进行插值
        float3 n0 = SimNormals[b.SimI0];
        float3 n1 = SimNormals[b.SimI1];
        float3 n2 = SimNormals[b.SimI2];

        float3 interpNormal = b.U * n0 + b.V * n1 + b.W * n2;
        float len = math.length(interpNormal);
        RenderNormals[i] = len > 1e-8f ? interpNormal / len : new float3(0, 1, 0);
    }
}

/// <summary>
/// 计算模拟网格顶点法线的Job（面积加权平均法线）
/// 这个法线是"平滑法线"——每个顶点的法线是其相邻面法线的加权平均
/// 用于渲染时产生平滑的光照效果
/// 
/// 三角形索引以 int3 形式传入，可直接从 SurfaceTriangleIndex / TriangleIndex Buffer
/// 零拷贝 reinterpret 得到（两者皆为 {int,int,int} 同布局）。
/// </summary>
[BurstCompile]
public struct ComputeSimNormalsJob : IJob
{
    [ReadOnly] public NativeArray<float3> Positions;
    [ReadOnly] public NativeArray<int3> Triangles;
    public NativeArray<float3> Normals;

    public void Execute()
    {
        // 清零
        for (int i = 0; i < Normals.Length; i++)
            Normals[i] = float3.zero;

        // 累加面积加权法线
        int numTris = Triangles.Length;
        for (int t = 0; t < numTris; t++)
        {
            int3 tri = Triangles[t];
            int i0 = tri.x;
            int i1 = tri.y;
            int i2 = tri.z;

            float3 p0 = Positions[i0];
            float3 p1 = Positions[i1];
            float3 p2 = Positions[i2];

            float3 faceNormal = math.cross(p1 - p0, p2 - p0);

            Normals[i0] = Normals[i0] + faceNormal;
            Normals[i1] = Normals[i1] + faceNormal;
            Normals[i2] = Normals[i2] + faceNormal;
        }

        // 归一化
        for (int i = 0; i < Normals.Length; i++)
        {
            float len = math.length(Normals[i]);
            if (len > 1e-8f)
                Normals[i] = Normals[i] / len;
            else
                Normals[i] = new float3(0, 1, 0);
        }
    }
}
