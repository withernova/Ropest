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
    public float EdgeStiffness;    // EdgeConstraint的stiff
    public float BendTwistKs;      // BendingAndTwistingConstraint的ks
    public float Damping;          // 速度阻尼系数
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
