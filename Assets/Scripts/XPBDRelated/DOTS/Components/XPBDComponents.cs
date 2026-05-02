using Unity.Entities;
using Unity.Mathematics;


[InternalBufferCapacity(0)]
public struct ParticlePosition : IBufferElementData
{
    public float3 Value;
}

[InternalBufferCapacity(0)]
public struct ParticlePrevPosition : IBufferElementData
{
    public float3 Value;
}

[InternalBufferCapacity(0)]
public struct ParticleVelocity : IBufferElementData
{
    public float3 Value;
}

[InternalBufferCapacity(0)]
public struct ParticleInvMass : IBufferElementData
{
    public float Value;
}


// 粒子是否为"表面粒子"：仅表面粒子参与跨体碰撞（CrossBodyCollisionSystem）。
// - 布料：所有粒子都是表面（Value=1）
// - 软体：仅真正处于表面的顶点（属于任一表面三角形）为 1，内部粒子为 0
// 不加该 Buffer 的实体（Rope 等）在跨体碰撞里视作全表面。
[InternalBufferCapacity(0)]
public struct ParticleSurfaceFlag : IBufferElementData
{
    public byte Value; // 0 = 内部 / 不参与跨体碰撞，1 = 表面
}


[InternalBufferCapacity(0)]
public struct GhostPosition : IBufferElementData
{
    public float3 Value;
}

[InternalBufferCapacity(0)]
public struct GhostPrevPosition : IBufferElementData
{
    public float3 Value;
}

[InternalBufferCapacity(0)]
public struct GhostVelocity : IBufferElementData
{
    public float3 Value;
}

[InternalBufferCapacity(0)]
public struct GhostInvMass : IBufferElementData
{
    public float Value;
}

[InternalBufferCapacity(0)]
public struct RopeRestLength : IBufferElementData
{
    public float Value;
}

[InternalBufferCapacity(0)]
public struct InitDarbouxVector : IBufferElementData
{
    public float3 Value;
}


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


[InternalBufferCapacity(0)]
public struct BendTwistLambda : IBufferElementData
{
    public float3 Value;
}


// ============================================================================
// XPBD 共享 Buffer 组件（布料与软体共用）
// 说明：原 ClothEdge / SoftBodyEdge 字段完全一致，合并为 XPBDEdge；
//       原 ClothDistanceLambda / SoftBodyDistanceLambda 同理合并为 XPBDDistanceLambda。
// ============================================================================
[InternalBufferCapacity(0)]
public struct XPBDEdge : IBufferElementData
{
    public int IndexA;
    public int IndexB;
    public float RestLength;
}

[InternalBufferCapacity(0)]
public struct XPBDDistanceLambda : IBufferElementData
{
    public float Value;
}

// ============================================================================
// 图着色（Graph Coloring）相关 Buffer
// 说明：当启用图着色并行方案时，Spawner 会先对 XPBDEdge / Tetrahedron 做
//       贪心染色并按颜色重排：同一颜色组内的约束互不共享粒子，
//       可以安全地用 IJobParallelFor 并行求解。
//       XPBDEdgeColorRange 的每个元素表示一个颜色组在 XPBDEdge Buffer 中的
//       [Start, Start+Count) 连续区间。
// ============================================================================
[InternalBufferCapacity(0)]
public struct XPBDEdgeColorRange : IBufferElementData
{
    public int Start;
    public int Count;
}

// 四面体体积约束的颜色分组区间（软体专属）
[InternalBufferCapacity(0)]
public struct SoftBodyTetColorRange : IBufferElementData
{
    public int Start;
    public int Count;
}

[InternalBufferCapacity(0)]
public struct TriangleIndex : IBufferElementData
{
    public int I0;
    public int I1;
    public int I2;
}


[InternalBufferCapacity(0)]
public struct SectionVertexIndex : IBufferElementData
{
    public int Value;
}


[InternalBufferCapacity(0)]
public struct MeshVertex : IBufferElementData
{
    public float3 Value;
}


public struct RopeSolverConfig : IComponentData
{
    public int NumPoints;
    public int NumGhostPoints;
    public int Segments;
    public int Subdivision;
    public int NumSubSteps;
    public float Radius;
    public float GhostDistance;
    public float GravityFactor;
    public float3 Gravity;
    public float EdgeStiffness;
    public float BendTwistStiffness;
    public float BendTwistKs;
    public float Damping;
    public float Friction;
}

// ============================================================================
// XPBD 基础共享配置（布料/软体共用字段）
// DOTS 的 IComponentData 为 struct 无法继承，这里通过"组合"让两种配置内嵌同一个
// XPBDBaseConfig，公共流水线 Job（PreSolve/PostSolve/DistanceConstraint/AnalyticalCollision）
// 只读取 Base 部分即可，专属字段（如 Subdivision / VolumeStiffness）各自保留。
// ============================================================================
public struct XPBDBaseConfig
{
    public int NumParticles;
    public int NumSubSteps;
    public float3 Gravity;
    public float DistanceStiffness;
    public float Damping;
    public float CollisionRadius;
    public float Friction;
}

public struct ClothSolverConfig : IComponentData
{
    public XPBDBaseConfig Base;
    public int Subdivision;
    public bool SkipAnalyticalCollision;
    public bool SkipCrossBodyCollision;

    // 是否使用基于图着色的并行距离约束求解。
    // 开启前 Spawner 必须对 XPBDEdge 按颜色重排并填充 XPBDEdgeColorRange Buffer。
    public bool UseGraphColoring;

    // === 向后兼容字段（旧代码 cfg.NumParticles / cfg.Gravity 等写法依然可用） ===
    public int NumParticles { get => Base.NumParticles; set => Base.NumParticles = value; }
    public int NumSubSteps { get => Base.NumSubSteps; set => Base.NumSubSteps = value; }
    public float3 Gravity { get => Base.Gravity; set => Base.Gravity = value; }
    public float DistanceStiffness { get => Base.DistanceStiffness; set => Base.DistanceStiffness = value; }
    public float Damping { get => Base.Damping; set => Base.Damping = value; }
    public float CollisionRadius { get => Base.CollisionRadius; set => Base.CollisionRadius = value; }
    public float Friction { get => Base.Friction; set => Base.Friction = value; }
}

public struct MeshUpdateTag : IComponentData { }
public struct RopeTag : IComponentData { }
public struct ClothTag : IComponentData { }

public enum AnalyticalColliderType
{
    Sphere = 0,
    Box = 1
}

public struct AnalyticalColliderData
{
    public AnalyticalColliderType Type;
    public float3 Center;
    public float Radius;
    public float3 HalfExtents;
    public quaternion Rotation;
    public quaternion InvRotation;
}

public struct SoftBodyTag : IComponentData { }

public struct SoftBodySolverConfig : IComponentData
{
    public XPBDBaseConfig Base;
    public float VolumeStiffness;

    // 是否使用基于图着色的并行约束求解（距离 + 体积都会启用）。
    // 开启前 Spawner 必须对 XPBDEdge / Tetrahedron 按颜色重排，
    // 并分别填充 XPBDEdgeColorRange / SoftBodyTetColorRange。
    public bool UseGraphColoring;

    // === 向后兼容字段 ===
    public int NumParticles { get => Base.NumParticles; set => Base.NumParticles = value; }
    public int NumSubSteps { get => Base.NumSubSteps; set => Base.NumSubSteps = value; }
    public float3 Gravity { get => Base.Gravity; set => Base.Gravity = value; }
    public float DistanceStiffness { get => Base.DistanceStiffness; set => Base.DistanceStiffness = value; }
    public float Damping { get => Base.Damping; set => Base.Damping = value; }
    public float CollisionRadius { get => Base.CollisionRadius; set => Base.CollisionRadius = value; }
    public float Friction { get => Base.Friction; set => Base.Friction = value; }
}

[InternalBufferCapacity(0)]
public struct Tetrahedron : IBufferElementData
{
    public int I0;
    public int I1;
    public int I2;
    public int I3;
}

[InternalBufferCapacity(0)]
public struct TetrahedronRestVolume : IBufferElementData
{
    public float Value;
}

[InternalBufferCapacity(0)]
public struct TetrahedronVolumeLambda : IBufferElementData
{
    public float Value;
}

[InternalBufferCapacity(0)]
public struct SurfaceTriangleIndex : IBufferElementData
{
    public int I0;
    public int I1;
    public int I2;
}

[InternalBufferCapacity(0)]
public struct RenderVertexBinding : IBufferElementData
{
    public int SimI0;
    public int SimI1;
    public int SimI2;
    public float U;
    public float V;
    public float W;
}

public struct RenderMeshConfig : IComponentData
{
    public int NumRenderVertices;
    public bool UseSubdivision;
}
