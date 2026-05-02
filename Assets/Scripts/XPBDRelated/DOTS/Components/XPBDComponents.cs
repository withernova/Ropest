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
