using System.Collections;
using System.Collections.Generic;
using System.IO;
using UnityEngine;
using UnityEngine.Profiling; // Recorder

/// <summary>
/// 综合应用场景 Benchmark（对应论文 5.4 节 BasketballGame）。
///
/// 工作流程：
/// - 直接挂在任意常驻 GameObject 上，启动后脚本会：
///   1. 等 <see cref="warmupSeconds"/> 秒（让场景 spawn 完成、玩家操作介入后进入稳态）；
///   2. 在随后的 <see cref="durationSeconds"/> 秒内每帧采集 Cloth / SoftBody 的
///      <c>XPBD.*.Frame</c> 与 <c>XPBD.*.Constraint</c> 两个 ProfilerMarker，
///      同时用 unscaledDeltaTime 记录端到端 FPS；
///   3. 结束后把每秒均值、整段 p50/p99、峰值粒子数（通过场景里
///      <see cref="CountActiveParticles"/> 反射计数）写入 CSV。
///
/// 输出两个 CSV：
/// - {outputCsvPath}          ：汇总统计（1 行）
/// - {outputCsvPath}.timeline ：每秒一行的时间序列，方便画 FPS/耗时曲线
/// </summary>
[DefaultExecutionOrder(-10000)]
public class BasketballGameBenchmarkRunner : MonoBehaviour
{
    [Header("采样时长（秒）")]
    [Min(1)] public float warmupSeconds = 5f;
    [Min(1)] public float durationSeconds = 60f;

    [Header("Marker 源")]
    [Tooltip("同时采集 Cloth 与 SoftBody 两组 Marker，缺席的那一组会是 0。")]
    public bool samplePhysicsMarkers = true;

    [Header("输出")]
    public string outputCsvPath = "xpbd_basketball_benchmark.csv";

    // ---------------- 运行时状态 ----------------
    private Recorder _clothFrame, _clothConstraint;
    private Recorder _softFrame, _softConstraint;

    IEnumerator Start()
    {
        // 综合场景保持正常游戏 VSync 条件也可；这里强制关闭 VSync 以测出真实负载上限
        QualitySettings.vSyncCount = 0;
        Application.targetFrameRate = -1;

        // 启用 Cloth/SoftBody 的同步采样（Marker 才能量到真实 Job 耗时）
        ClothSimulationSystem.BenchmarkSyncMode = true;
        SoftBodySimulationSystem.BenchmarkSyncMode = true;

        if (samplePhysicsMarkers)
        {
            _clothFrame      = Recorder.Get("XPBD.Cloth.Frame");
            _clothConstraint = Recorder.Get("XPBD.Cloth.Constraint");
            _softFrame       = Recorder.Get("XPBD.SoftBody.Frame");
            _softConstraint  = Recorder.Get("XPBD.SoftBody.Constraint");
            _clothFrame.enabled = _clothConstraint.enabled = true;
            _softFrame.enabled  = _softConstraint.enabled  = true;
        }

        Debug.Log($"[BBallBench] 预热 {warmupSeconds}s ...");
        float t0 = Time.unscaledTime;
        while (Time.unscaledTime - t0 < warmupSeconds) yield return null;

        Debug.Log($"[BBallBench] 开始采集 {durationSeconds}s ...");

        var dtAll       = new List<double>(4096);
        var clothFAll   = new List<double>(4096);
        var clothCAll   = new List<double>(4096);
        var softFAll    = new List<double>(4096);
        var softCAll    = new List<double>(4096);
        var timeline    = new List<string>(120);

        // 每秒桶累计
        double bucketTs = Time.unscaledTime;
        int bucketFrames = 0;
        double bucketClothF = 0, bucketClothC = 0, bucketSoftF = 0, bucketSoftC = 0, bucketDt = 0;
        int peakParticles = 0;

        timeline.Add("sec,fps_avg,frame_avg_ms,cloth_frame_ms,cloth_cons_ms,soft_frame_ms,soft_cons_ms,active_particles");

        float tStart = Time.unscaledTime;
        while (Time.unscaledTime - tStart < durationSeconds)
        {
            yield return null;

            double dt = Time.unscaledDeltaTime;
            double cF = _clothFrame      != null ? _clothFrame.elapsedNanoseconds      * 1e-6 : 0;
            double cC = _clothConstraint != null ? _clothConstraint.elapsedNanoseconds * 1e-6 : 0;
            double sF = _softFrame       != null ? _softFrame.elapsedNanoseconds       * 1e-6 : 0;
            double sC = _softConstraint  != null ? _softConstraint.elapsedNanoseconds  * 1e-6 : 0;

            if (dt > 0) dtAll.Add(dt);
            if (cF > 0) clothFAll.Add(cF);
            if (cC > 0) clothCAll.Add(cC);
            if (sF > 0) softFAll.Add(sF);
            if (sC > 0) softCAll.Add(sC);

            bucketFrames++;
            bucketDt      += dt;
            bucketClothF  += cF;
            bucketClothC  += cC;
            bucketSoftF   += sF;
            bucketSoftC   += sC;

            // 每 1 秒统计一次 & 写一行时间序列
            if (Time.unscaledTime - bucketTs >= 1.0)
            {
                int activeParticles = CountActiveParticles();
                if (activeParticles > peakParticles) peakParticles = activeParticles;

                double avgDt = bucketFrames > 0 ? bucketDt / bucketFrames : 0;
                double avgFps = avgDt > 0 ? 1.0 / avgDt : 0;

                int secIdx = timeline.Count; // 从 1 开始正好，因为第 0 行是表头
                timeline.Add(string.Format(
                    "{0},{1:F1},{2:F2},{3:F2},{4:F2},{5:F2},{6:F2},{7}",
                    secIdx,
                    avgFps,
                    avgDt * 1000,
                    bucketFrames > 0 ? bucketClothF / bucketFrames : 0,
                    bucketFrames > 0 ? bucketClothC / bucketFrames : 0,
                    bucketFrames > 0 ? bucketSoftF  / bucketFrames : 0,
                    bucketFrames > 0 ? bucketSoftC  / bucketFrames : 0,
                    activeParticles));

                bucketTs = Time.unscaledTime;
                bucketFrames = 0;
                bucketDt = bucketClothF = bucketClothC = bucketSoftF = bucketSoftC = 0;
            }
        }

        // 汇总
        var (dtMean, dtP50, dtP99) = Stats(dtAll);
        var (cFm, cFp50, cFp99)    = Stats(clothFAll);
        var (cCm, cCp50, cCp99)    = Stats(clothCAll);
        var (sFm, sFp50, sFp99)    = Stats(softFAll);
        var (sCm, sCp50, sCp99)    = Stats(softCAll);
        double fpsMean = dtP50 > 0 ? 1.0 / dtP50 : 0;
        double fpsP99Low = dtP99 > 0 ? 1.0 / dtP99 : 0; // 1% Low 概念：p99 帧时间倒数

        var summary = new System.Text.StringBuilder();
        summary.AppendLine("metric,value");
        summary.AppendLine($"duration_s,{durationSeconds:F1}");
        summary.AppendLine($"sample_frames,{dtAll.Count}");
        summary.AppendLine($"fps_mean,{fpsMean:F1}");
        summary.AppendLine($"fps_1pct_low,{fpsP99Low:F1}");
        summary.AppendLine($"frame_ms_mean,{dtMean*1000:F2}");
        summary.AppendLine($"frame_ms_p50,{dtP50*1000:F2}");
        summary.AppendLine($"frame_ms_p99,{dtP99*1000:F2}");
        summary.AppendLine($"cloth_frame_ms_p50,{cFp50:F2}");
        summary.AppendLine($"cloth_cons_ms_p50,{cCp50:F2}");
        summary.AppendLine($"soft_frame_ms_p50,{sFp50:F2}");
        summary.AppendLine($"soft_cons_ms_p50,{sCp50:F2}");
        summary.AppendLine($"peak_active_particles,{peakParticles}");

        var summaryPath = Path.Combine(Application.persistentDataPath, outputCsvPath);
        File.WriteAllText(summaryPath, summary.ToString());

        var timelinePath = summaryPath + ".timeline.csv";
        File.WriteAllText(timelinePath, string.Join("\n", timeline));

        Debug.Log($"[BBallBench] 汇总: FPS {fpsMean:F1} (1% low {fpsP99Low:F1}), " +
                  $"peak particles {peakParticles}, cloth {cFp50:F2}/{cCp50:F2} ms, " +
                  $"soft {sFp50:F2}/{sCp50:F2} ms");
        Debug.Log($"[BBallBench] summary → {summaryPath}");
        Debug.Log($"[BBallBench] timeline → {timelinePath}");
    }

    /// <summary>
    /// 扫描场景里所有 Cloth DOTS 实体与 SoftBody DOTS 实体的粒子总数。
    /// 通过遍历 World 的 ParticlePosition Buffer 完成（无反射，Burst 编译无影响）。
    /// </summary>
    private static int CountActiveParticles()
    {
        int total = 0;
        var world = Unity.Entities.World.DefaultGameObjectInjectionWorld;
        if (world == null || !world.IsCreated) return 0;

        // 布料粒子
        var em = world.EntityManager;
        using (var q = em.CreateEntityQuery(
                   Unity.Entities.ComponentType.ReadOnly<ClothTag>(),
                   Unity.Entities.ComponentType.ReadOnly<ParticlePosition>()))
        {
            using (var entities = q.ToEntityArray(Unity.Collections.Allocator.Temp))
            {
                foreach (var e in entities)
                    total += em.GetBuffer<ParticlePosition>(e).Length;
            }
        }

        // 软体粒子
        using (var q = em.CreateEntityQuery(
                   Unity.Entities.ComponentType.ReadOnly<SoftBodyTag>(),
                   Unity.Entities.ComponentType.ReadOnly<ParticlePosition>()))
        {
            using (var entities = q.ToEntityArray(Unity.Collections.Allocator.Temp))
            {
                foreach (var e in entities)
                    total += em.GetBuffer<ParticlePosition>(e).Length;
            }
        }

        return total;
    }

    static (double mean, double p50, double p99) Stats(List<double> xs)
    {
        if (xs == null || xs.Count == 0) return (0, 0, 0);
        xs.Sort();
        double sum = 0; foreach (var v in xs) sum += v;
        double mean = sum / xs.Count;
        double p50 = xs[xs.Count / 2];
        int i99 = Mathf.Clamp(Mathf.FloorToInt(xs.Count * 0.99f), 0, xs.Count - 1);
        double p99 = xs[i99];
        return (mean, p50, p99);
    }
}
