using Unity.Entities;
using UnityEngine;

/// <summary>
/// XPBD 物理仿真的全局调试控制器。
/// - IsPaused：当前是否处于暂停状态。暂停后 ClothSimulationSystem / SoftBodySimulationSystem 
///   将跳过每帧的 PreSolve → SubSteps → PostSolve 流程，渲染系统继续运行以保持画面可观察。
/// - PendingSubSteps：暂停状态下外部请求推进的"子步迭代次数"。一次 FixedStep 消费 1 个计数，
///   执行一次"NumSubSteps=1"的完整帧（PreSolve + 1 个约束迭代 + PostSolve），
///   从而让用户可以按键逐迭代地查看约束求解中间结果。
/// 
/// 为什么用静态类：DOTS 的 ISystem 是 struct，不适合承载托管状态；静态类对 MonoBehaviour/ECS 
/// 双侧访问最简单，且本身就是单例调试工具，无需支持多实例。
/// </summary>
public static class XPBDDebugController
{
    /// <summary>是否暂停物理仿真（仅影响布料和软体的 XPBD 主循环）。</summary>
    public static bool IsPaused { get; private set; }

    /// <summary>待执行的"逐迭代推进"请求计数。</summary>
    public static int PendingSubSteps { get; private set; }

    /// <summary>
    /// 当前 FixedStep 是否被允许执行物理仿真。
    /// 由 XPBDDebugTickSystem 在每次 FixedStep 开头写入（仅一次），
    /// Cloth / SoftBody 模拟系统读它来决定是否跳过 OnUpdate。
    /// 这样同一 FixedStep 内的多个系统看到同一个值，计数不会被重复消耗。
    /// </summary>
    public static bool AllowCurrentFixedStep { get; internal set; } = true;

    /// <summary>是否处于"单步模式"：本 FixedStep 是因 PendingSubSteps 推进而被允许的。</summary>
    public static bool IsSingleStepFrame { get; internal set; }

    /// <summary>
    /// 从进入游戏到现在、物理系统总共实际推进过的 FixedStep 数量（正常 + 单步都会累加）。
    /// 跳过的（暂停且无步进请求）不计入。这是一个"全局逻辑帧号"，方便定位。
    /// </summary>
    public static int TotalPhysicsTicks { get; private set; }

    /// <summary>
    /// 当前暂停会话中，已通过 N/M 键推进的子步序号。
    /// 从 1 开始计数；每次"恢复正常运行（P 取消暂停）"时归零。
    /// 这是最直观的"当前是第几个迭代"——看这个就知道从按 P 暂停到现在一共按了多少次推进。
    /// </summary>
    public static int SingleStepIndex { get; private set; }

    /// <summary>
    /// 最近一次 FixedStep 实际执行的 SubStep 次数。
    /// 正常运行 = 配置的 NumSubSteps；单步模式 = 1；被跳过的帧不更新。
    /// </summary>
    public static int LastSubStepsExecuted { get; internal set; }

    /// <summary>是否开启 Gizmos 速度向量绘制。</summary>
    public static bool DrawVelocityGizmos { get; set; } = true;

    /// <summary>Gizmos 速度向量的缩放（米/秒 → 场景米）。</summary>
    public static float VelocityGizmoScale { get; set; } = 0.1f;

    public static void SetPaused(bool paused)
    {
        IsPaused = paused;
        if (!paused)
        {
            // 恢复正常运行：清掉未消费的步进请求 + 复位单步序号（进入下一段暂停会话时从 1 重新计数）。
            PendingSubSteps = 0;
            SingleStepIndex = 0;
        }
    }

    public static void TogglePause() => SetPaused(!IsPaused);

    /// <summary>外部（HUD / 输入脚本）请求推进 n 次子步。</summary>
    public static void RequestStep(int count = 1)
    {
        if (!IsPaused) return;
        if (count <= 0) return;
        PendingSubSteps += count;
    }

    /// <summary>
    /// 内部：FixedStep 开始时由 XPBDDebugTickSystem 调用。
    /// 返回本次 FixedStep 是否应当执行物理。
    /// 同时在这里维护所有"迭代计数器"——由于此方法在 FixedStep 最前面只会被调用一次，
    /// 累加不会重复。
    /// </summary>
    internal static bool ConsumeFixedStepTick()
    {
        if (!IsPaused)
        {
            AllowCurrentFixedStep = true;
            IsSingleStepFrame = false;
            TotalPhysicsTicks++;
            // 正常运行时 SubStepsExecuted 由 Cloth/SoftBody 系统各自的配置决定，
            // 这里写一个"占位值 -1"告诉 UI "看 solver 配置"。
            LastSubStepsExecuted = -1;
            return true;
        }
        if (PendingSubSteps > 0)
        {
            PendingSubSteps--;
            AllowCurrentFixedStep = true;
            IsSingleStepFrame = true;
            TotalPhysicsTicks++;
            SingleStepIndex++;
            LastSubStepsExecuted = 1;
            return true;
        }
        AllowCurrentFixedStep = false;
        IsSingleStepFrame = false;
        return false;
    }
}

/// <summary>
/// 在 FixedStepSimulationSystemGroup 最前面执行，统一决定"本次 FixedStep 是否跑物理"。
/// OrderFirst=true 确保它早于 ClothSimulationSystem / SoftBodySimulationSystem。
/// </summary>
[UpdateInGroup(typeof(FixedStepSimulationSystemGroup), OrderFirst = true)]
public partial struct XPBDDebugTickSystem : ISystem
{
    public void OnUpdate(ref SystemState state)
    {
        XPBDDebugController.ConsumeFixedStepTick();
    }
}