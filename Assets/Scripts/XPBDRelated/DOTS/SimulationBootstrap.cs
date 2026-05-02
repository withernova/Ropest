using UnityEngine;

public static class SimulationBootstrap
{
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
