using Unity.Entities;
using Unity.Mathematics;

// ============================================================
// 粒子数据 Buffers
// ============================================================

/// <summary>
/// 粒子位置（Rope用pointPos，Cloth用pos）
/// </summary>
[InternalBufferCapacity(0)]
public struct ParticlePosition : IBufferElementData
{
    public float3 Value;
}

/// <summary>
/// 粒子前一帧位置
/// </summary>
[InternalBufferCapacity(0)]
public struct ParticlePrevPosition : IBufferElementData
{
    public float3 Value;
}

/// <summary>
/// 粒子速度
/// </summary>
[InternalBufferCapacity(0)]
public struct ParticleVelocity : IBufferElementData
{
    public float3 Value;
}

/// <summary>
/// 粒子逆质量
/// </summary>
[InternalBufferCapacity(0)]
public struct ParticleInvMass : IBufferElementData
{
    public float Value;
}

// ============================================================
// Rope 专用 Buffers
// ============================================================

/// <summary>
/// Rope的Ghost点位置
/// </summary>
[InternalBufferCapacity(0)]
public struct GhostPosition : IBufferElementData
{
    public float3 Value;
}

/// <summary>
/// Rope的Ghost点前一帧位置
/// </summary>
[InternalBufferCapacity(0)]
public struct GhostPrevPosition : IBufferElementData
{
    public float3 Value;
}

/// <summary>
/// Rope的Ghost点速度
/// </summary>
[InternalBufferCapacity(0)]
public struct GhostVelocity : IBufferElementData
{
    public float3 Value;
}

/// <summary>
/// Rope的Ghost点逆质量
/// </summary>
[InternalBufferCapacity(0)]
public struct GhostInvMass : IBufferElementData
{
    public float Value;
}

/// <summary>
/// Rope各段的静止长度
/// </summary>
[InternalBufferCapacity(0)]
public struct RopeRestLength : IBufferElementData
{
    public float Value;
}

/// <summary>
/// Rope的初始Darboux向量（BendingAndTwisting约束用）
/// </summary>
[InternalBufferCapacity(0)]
public struct InitDarbouxVector : IBufferElementData
{
    public float3 Value;
}

// ============================================================
// Rope EdgeConstraint Lambda Buffers
// ============================================================

[InternalBufferCapacity(0)]
public struct EdgeLambda0 : IBufferElementData
{
    public float Value;
}

[InternalBufferCapacity(0)]
public struct EdgeLambda1 : IBufferElementData
{
    public float Value;
}

[InternalBufferCapacity(0)]
public struct EdgeLambda2 : IBufferElementData
{
    public float Value;
}

// ============================================================
// Rope BendTwist Lambda Buffer
// ============================================================

[InternalBufferCapacity(0)]
public struct BendTwistLambda : IBufferElementData
{
    public float3 Value;
}

// ============================================================
// Cloth 专用 Buffers
// ============================================================

/// <summary>
/// Cloth距离约束的边数据
/// </summary>
[InternalBufferCapacity(0)]
public struct ClothEdge : IBufferElementData
{
    public int IndexA;
    public int IndexB;
    public float RestLength;
}

/// <summary>
/// Cloth距离约束的Lambda
/// </summary>
[InternalBufferCapacity(0)]
public struct ClothDistanceLambda : IBufferElementData
{
    public float Value;
}

/// <summary>
/// 三角形索引数据
/// </summary>
[InternalBufferCapacity(0)]
public struct TriangleIndex : IBufferElementData
{
    public int I0;
    public int I1;
    public int I2;
}

// ============================================================
// Rope 截面映射 Buffer（用于渲染：截面索引 -> 顶点索引）
// ============================================================

/// <summary>
/// 截面中的顶点索引（扁平化存储，每subdivision个为一组）
/// </summary>
[InternalBufferCapacity(0)]
public struct SectionVertexIndex : IBufferElementData
{
    public int Value;
}

// ============================================================
// Mesh顶点Buffer（用于渲染回写）
// ============================================================

[InternalBufferCapacity(0)]
public struct MeshVertex : IBufferElementData
{
    public float3 Value;
}

// ============================================================
// 配置组件
// ============================================================

/// <summary>
/// Rope求解器配置
/// </summary>
public struct RopeSolverConfig : IComponentData
{
    public int NumPoints;          // pointPos数量
    public int NumGhostPoints;     // ghostPos数量
    public int Segments;
    public int Subdivision;
    public int NumSubSteps;
    public float Radius;
    public float GhostDistance;
    public float GravityFactor;
    public float3 Gravity;
    public float EdgeStiffness;        // EdgeConstraint的compliance (α)，原版默认0.0001
    public float BendTwistStiffness;   // BendingAndTwistingConstraint的compliance (α)，原版默认0.6
    public float BendTwistKs;          // BendingAndTwistingConstraint的ks（刚度缩放系数），原版默认0.8
    public float Damping;              // 速度阻尼系数
    public float Friction;             // 碰撞摩擦系数（0~1），0=无摩擦，推荐 0.2~0.5
}

/// <summary>
/// Cloth求解器配置
/// </summary>
public struct ClothSolverConfig : IComponentData
{
    public int NumParticles;
    public int Subdivision;
    public int NumSubSteps;
    public float3 Gravity;
    public float DistanceStiffness;
    public float Damping;           // 速度阻尼系数（0~1，越大阻尼越强，推荐0.01~0.05）
    public float CollisionRadius;   // 碰撞检测半径（安全边距）
    public float Friction;          // 碰撞摩擦系数（0~1，推荐0.3~0.6）
    // 运行时开关：设为 true 时本帧跳过与场景解析碰撞体（地板/Box/球体）的碰撞
    // 用于"布料升起/落下"等非交互阶段，避免布料穿过地面产生被压扁的怪异网格
    public bool SkipAnalyticalCollision;
    // 运行时开关：设为 true 时本帧跳过与其他 Body（软体球等）的跨体碰撞
    public bool SkipCrossBodyCollision;
}

/// <summary>
/// 标记需要更新Mesh的Entity
/// </summary>
public struct MeshUpdateTag : IComponentData { }

/// <summary>
/// 标记Rope Entity
/// </summary>
public struct RopeTag : IComponentData { }

/// <summary>
/// 标记Cloth Entity
/// </summary>
public struct ClothTag : IComponentData { }

// ============================================================
// 解析碰撞体数据（用于替代Unity Physics的碰撞检测）
// 支持球体和Box，可Burst化
// ============================================================

/// <summary>
/// 碰撞体类型枚举
/// </summary>
public enum AnalyticalColliderType
{
    Sphere = 0,
    Box = 1
}

/// <summary>
/// 解析碰撞体数据（blittable，可在NativeArray中使用）
/// 球体：Center + Radius
/// Box：Center + HalfExtents + InverseRotation（用于将点变换到Box局部空间）
/// </summary>
public struct AnalyticalColliderData
{
    public AnalyticalColliderType Type;
    public float3 Center;
    public float Radius;            // 球体半径
    public float3 HalfExtents;      // Box半尺寸
    public quaternion Rotation;     // Box旋转
    public quaternion InvRotation;  // Box逆旋转（用于将世界坐标变换到局部坐标）
}

// ============================================================
// SoftBody（可变形软体）专用 Buffers 和 Components
// ============================================================

/// <summary>
/// 标记SoftBody Entity
/// </summary>
public struct SoftBodyTag : IComponentData { }

/// <summary>
/// SoftBody求解器配置
/// </summary>
public struct SoftBodySolverConfig : IComponentData
{
    public int NumParticles;        // 粒子（顶点）数量
    public int NumSubSteps;         // 子步迭代次数
    public float3 Gravity;          // 重力
    public float DistanceStiffness; // 距离约束柔度（compliance），0=刚性
    public float VolumeStiffness;   // 体积约束柔度（compliance），0=不可压缩
    public float Damping;           // 速度阻尼系数（0~1）
    public float CollisionRadius;   // 碰撞检测半径
    public float Friction;          // 碰撞摩擦系数
}

/// <summary>
/// SoftBody距离约束的边数据
/// </summary>
[InternalBufferCapacity(0)]
public struct SoftBodyEdge : IBufferElementData
{
    public int IndexA;
    public int IndexB;
    public float RestLength;
}

/// <summary>
/// SoftBody距离约束的Lambda
/// </summary>
[InternalBufferCapacity(0)]
public struct SoftBodyDistanceLambda : IBufferElementData
{
    public float Value;
}

/// <summary>
/// 四面体索引数据（4个顶点索引）
/// </summary>
[InternalBufferCapacity(0)]
public struct Tetrahedron : IBufferElementData
{
    public int I0;
    public int I1;
    public int I2;
    public int I3;
}

/// <summary>
/// 四面体静止体积
/// </summary>
[InternalBufferCapacity(0)]
public struct TetrahedronRestVolume : IBufferElementData
{
    public float Value;
}

/// <summary>
/// 四面体体积约束的Lambda
/// </summary>
[InternalBufferCapacity(0)]
public struct TetrahedronVolumeLambda : IBufferElementData
{
    public float Value;
}

/// <summary>
/// 表面三角形索引（用于渲染）
/// </summary>
[InternalBufferCapacity(0)]
public struct SurfaceTriangleIndex : IBufferElementData
{
    public int I0;
    public int I1;
    public int I2;
}

// ============================================================
// 渲染网格绑定数据（模拟网格 -> 高面数渲染网格的重心坐标插值）
// ============================================================

/// <summary>
/// 渲染顶点的重心坐标绑定数据
/// 每个渲染顶点绑定到一个模拟三角形上，通过重心坐标插值得到位置
/// </summary>
[InternalBufferCapacity(0)]
public struct RenderVertexBinding : IBufferElementData
{
    /// <summary>模拟三角形的三个顶点索引（对应ParticlePosition中的索引）</summary>
    public int SimI0;
    public int SimI1;
    public int SimI2;

    /// <summary>重心坐标权重 (u, v, w)，满足 u + v + w = 1</summary>
    public float U;
    public float V;
    public float W;
}

/// <summary>
/// 渲染网格配置（存储渲染顶点数量等信息）
/// </summary>
public struct RenderMeshConfig : IComponentData
{
    /// <summary>渲染网格顶点数量</summary>
    public int NumRenderVertices;
    /// <summary>是否启用高面数渲染（false则直接用模拟网格渲染）</summary>
    public bool UseSubdivision;
}
