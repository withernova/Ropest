using UnityEngine;

/// <summary>
/// 仿真相关的全局启动参数设置。
///
/// 【核心作用】限制 Unity 在卡顿后的 FixedUpdate 追帧雪崩：
///   Unity 默认 Time.maximumDeltaTime = 1/3 秒，意味着单帧最长可以追 16 个 fixed tick
///   （1/0.02 = 16 次 SimulationSystemGroup 执行）。当 SoftBody/Cloth/CrossBodyCollision
///   单次 tick 约 5ms 时，一旦本帧被卡到 20ms，下一帧就要补跑 3~4 次 tick，
///   直接把帧耗时放大到 40~50ms，形成卡顿-补帧-更卡顿的恶性循环。
///
/// 把上限压到 1/30s = 0.0333s，意味着最多只追 ~1.67 次 tick，主帧时间被硬性封顶，
/// 以"仿真时间略微变慢"为代价换取"帧率稳定不雪崩"。
///
/// 注意：这只改变 Unity 追赶真实时间的行为，不改变任何仿真逻辑/物理参数。
/// </summary>
public static class SimulationBootstrap
{
    // 一帧最多允许追赶的真实时间（秒）。
    // 0.0333f => 最多追约 1.67 个 fixed tick（fixedDeltaTime 默认 0.02）
    // 如果仍感到卡顿波动大，可再降到 0.025f（最多 1.25 次）。
    private const float MaxDeltaTime = 1f / 30f;

    [RuntimeInitializeOnLoadMethod(RuntimeInitializeLoadType.BeforeSceneLoad)]
    private static void ApplyTimeSettings()
    {
        Time.maximumDeltaTime = MaxDeltaTime;
        Debug.Log($"[SimulationBootstrap] Time.maximumDeltaTime = {Time.maximumDeltaTime:F4}s " +
                  $"(fixedDeltaTime = {Time.fixedDeltaTime:F4}s, 每帧最多追 " +
                  $"{Time.maximumDeltaTime / Time.fixedDeltaTime:F2} 个 fixed tick)");
    }
}
