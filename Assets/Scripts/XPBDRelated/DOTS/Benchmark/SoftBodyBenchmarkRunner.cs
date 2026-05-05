using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Text;
using Unity.Jobs.LowLevel.Unsafe;
using UnityEngine;
using UnityEngine.Profiling; // Recorder

/// <summary>
/// 软体仿真 Benchmark 自动采数脚本（对应论文 5.3 节）。
///
/// 与布料 <see cref="BenchmarkRunner"/> 的差异：
/// - 驱动对象为 <see cref="SoftBodyRuntimeSpawner"/>；
/// - 扫描维度为 软体体素分辨率 <c>resolution</c>（对应内部粒子/四面体密度）；
/// - ProfilerMarker 为 XPBD.SoftBody.Frame / XPBD.SoftBody.Constraint；
/// - 同样输出 CSV，用法一致。
///
/// 运行前确认：
/// 1. 场景内有一个 <see cref="SoftBodyRuntimeSpawner"/>，把它拖进 targetSpawner；
/// 2. <see cref="SoftBodySimulationSystem.BenchmarkSyncMode"/> 会被脚本自动开启。
/// </summary>
[DefaultExecutionOrder(-10000)]
public class SoftBodyBenchmarkRunner : MonoBehaviour
{
    // ---------------- Inspector 配置 ----------------

    [Header("被测 Spawner（模板）")]
    [Tooltip("场景中的 SoftBodyRuntimeSpawner 引用。会被 Destroy 并按配置重新 Instantiate。")]
    public SoftBodyRuntimeSpawner targetSpawner;

    [Header("采样策略")]
    [Min(60)] public int framesPerCase = 600;
    [Min(0)] public int warmupFrames = 120;
    [Min(0)] public int coolDownFrames = 30;

    [Header("实验 A：软体分辨率扫描")]
    public bool runResolutionSweep = true;
    [Tooltip("体素分辨率列表。体素粒子数 ≈ (r+1)³，经球面投影后有效粒子约 52%%。\n" +
             "默认 6/10/16/24/36/50 → 体素粒子 343/1331/4.9k/15.6k/50.7k/132.7k，\n" +
             "量级跨度与布料 Benchmark（32²~256² 即 1k~65k 粒子）相当。\n" +
             "注意 r≥36 会产生 5×r³ ≈ 18万+ 四面体约束，显存/内存压力显著，测试前先确认机器容量。")]
    public int[] resolutionCases = new[] { 6, 10, 16, 24, 36, 50 };

    [Header("实验 B：图着色 vs 串行 IJob")]
    public bool runColoringCompare = true;
    [Tooltip("图着色对比固定的分辨率。建议中等规模（24 ≈ 1.5万粒子）避开小规模调度开销噪声。")]
    public int coloringFixedResolution = 24;

    [Header("实验 C：多核扩展性")]
    public bool runThreadScaling = true;
    public int[] threadCases = new[] { 1, 2, 4, 6, 8, 12, 16 };
    [Tooltip("多核扫描固定分辨率。24 对应约 1.5 万粒子，足以压满多线程。")]
    public int threadFixedResolution = 24;

    [Header("输出")]
    [Tooltip("相对 Application.persistentDataPath 的 CSV 文件名。")]
    public string outputCsvPath = "xpbd_softbody_benchmark.csv";

    // ---------------- 运行时状态 ----------------

    private readonly StringBuilder _csv = new StringBuilder();
    private Recorder _frameRecorder;
    private Recorder _constraintRecorder;

    private GameObject _templateGO;
    private int _tplResolution;
    private bool _tplUseColoring;

    private SoftBodyRuntimeSpawner _activeSpawner;

    IEnumerator Start()
    {
        if (targetSpawner == null)
        {
            Debug.LogError("[SoftBench] targetSpawner 未设置，无法启动。");
            yield break;
        }

        // 模板冻结：同布料 Runner 思路
        _templateGO = targetSpawner.gameObject;
        _tplResolution = targetSpawner.resolution;
        _tplUseColoring = targetSpawner.useGraphColoring;
        _templateGO.SetActive(false);
        _activeSpawner = null;
        targetSpawner = null;

        // 解除帧率限制
        QualitySettings.vSyncCount = 0;
        Application.targetFrameRate = -1;

        // 打开软体 System 的同步采样模式
        SoftBodySimulationSystem.BenchmarkSyncMode = true;

        // CSV 表头（与布料 Runner 同构，方便后续统一后处理）
        _csv.AppendLine("experiment,case,resolution,numParticles,useColoring,threads," +
                        "frame_ms_mean,frame_ms_p50,frame_ms_p99," +
                        "constraint_ms_mean,constraint_ms_p50,constraint_ms_p99," +
                        "fps_mean");

        _frameRecorder = Recorder.Get("XPBD.SoftBody.Frame");
        _constraintRecorder = Recorder.Get("XPBD.SoftBody.Constraint");
        _frameRecorder.enabled = true;
        _constraintRecorder.enabled = true;

        // 预热
        yield return SpawnAndWait(_tplResolution, _tplUseColoring, JobsUtility.JobWorkerCount,
                                  warmupFrames + 30, collect: false, experiment: "warmup");

        // === 实验 A：分辨率扫描（每档跑 关/开 着色 两轮）===
        if (runResolutionSweep)
        {
            foreach (var useColor in new[] { false, true })
            {
                foreach (var r in resolutionCases)
                {
                    yield return SpawnAndWait(r, useColor,
                                              JobsUtility.JobWorkerCount,
                                              framesPerCase, collect: true,
                                              experiment: "resolution");
                    yield return CoolDown();
                }
            }
        }

        // === 实验 B：图着色对比 ===
        if (runColoringCompare)
        {
            foreach (var useColor in new[] { false, true })
            {
                yield return SpawnAndWait(coloringFixedResolution, useColor,
                                          JobsUtility.JobWorkerCount,
                                          framesPerCase, collect: true,
                                          experiment: "coloring");
                yield return CoolDown();
            }
        }

        // === 实验 C：多核扩展性 ===
        if (runThreadScaling)
        {
            int maxWorkers = JobsUtility.JobWorkerMaximumCount;
            int originalWorkers = JobsUtility.JobWorkerCount;
            foreach (var t in threadCases)
            {
                int clamped = Mathf.Clamp(t, 1, maxWorkers);
                JobsUtility.JobWorkerCount = clamped;
                yield return SpawnAndWait(threadFixedResolution, true, clamped,
                                          framesPerCase, collect: true,
                                          experiment: "threads");
                yield return CoolDown();
            }
            JobsUtility.JobWorkerCount = originalWorkers;
        }

        var path = Path.Combine(Application.persistentDataPath, outputCsvPath);
        File.WriteAllText(path, _csv.ToString());
        Debug.Log($"[SoftBench] CSV 已写入: {path}");
    }

    IEnumerator SpawnAndWait(int resolution, bool useColoring, int threadsInUse,
                              int frameCount, bool collect, string experiment)
    {
        if (_activeSpawner != null)
        {
            Destroy(_activeSpawner.gameObject);
            _activeSpawner = null;
            yield return null;
            yield return null;
        }

        if (_templateGO == null)
        {
            Debug.LogError("[SoftBench] 模板 GameObject 已丢失，无法继续。");
            yield break;
        }

        var go = Instantiate(_templateGO);
        go.name = $"SoftBench_res{resolution}_color{useColoring}";
        var sp = go.GetComponent<SoftBodyRuntimeSpawner>();
        sp.resolution = Mathf.Max(1, resolution);
        sp.useGraphColoring = useColoring;
        sp.enabled = true;
        go.SetActive(true);
        _activeSpawner = sp;

        // 等待 Spawn 完成（软体 mesh 构建比布料重，多等几帧）
        yield return null;
        yield return null;
        yield return null;
        yield return null;

        int warm = Mathf.Max(warmupFrames, 30);
        for (int i = 0; i < warm; i++) yield return null;

        if (!collect) yield break;

        var frameSamples = new List<double>(frameCount);
        var constraintSamples = new List<double>(frameCount);
        var frameDtSamples = new List<double>(frameCount);
        for (int i = 0; i < frameCount; i++)
        {
            yield return null;
            double frameMs = _frameRecorder.elapsedNanoseconds * 1e-6;
            double consMs = _constraintRecorder.elapsedNanoseconds * 1e-6;
            if (frameMs > 0) frameSamples.Add(frameMs);
            if (consMs > 0) constraintSamples.Add(consMs);
            if (Time.unscaledDeltaTime > 0) frameDtSamples.Add(Time.unscaledDeltaTime);
        }

        // 粒子数估算：体素立方体 (r+1)^3（实际会被球面投影过滤，但量级一致，用于 log/图表 x 轴）
        int numParticlesApprox = (resolution + 1) * (resolution + 1) * (resolution + 1);

        var (fm, fp50, fp99) = Stats(frameSamples);
        var (cm, cp50, cp99) = Stats(constraintSamples);
        var (_, dtMedian, _) = Stats(frameDtSamples);
        float fpsMean = dtMedian > 0 ? (float)(1.0 / dtMedian) : 0f;

        _csv.AppendLine($"{experiment},res{resolution},{resolution},{numParticlesApprox}," +
                        $"{useColoring},{threadsInUse}," +
                        $"{fm:F3},{fp50:F3},{fp99:F3},{cm:F3},{cp50:F3},{cp99:F3},{fpsMean:F1}");

        Debug.Log($"[SoftBench] {experiment} res={resolution} color={useColoring} " +
                  $"threads={threadsInUse} | frame {fm:F2}/{fp50:F2}/{fp99:F2} ms | " +
                  $"cons {cm:F2}/{cp50:F2}/{cp99:F2} ms | fps={fpsMean:F1}");
    }

    IEnumerator CoolDown()
    {
        for (int i = 0; i < coolDownFrames; i++) yield return null;
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
