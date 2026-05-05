using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Text;
using Unity.Jobs.LowLevel.Unsafe;
using Unity.Profiling;
using UnityEngine;
using UnityEngine.Profiling; // Recorder 所在命名空间

/// <summary>
/// XPBD DOTS 实验自动采数脚本（论文第 5 章用）。
///
/// 工作流程：
/// 1. 在 Inspector 把当前场景里的 <see cref="ClothRuntimeSpawner"/> 拖进 <see cref="targetSpawner"/>；
///    该 Spawner 所在的 GameObject 会被当作"模板"反复销毁/重建以切换分辨率。
/// 2. 在 Inspector 勾选要跑的实验项（分辨率扫描 / 图着色对比 / 多核扫描）。
/// 3. 运行场景，脚本会自动遍历所有配置，每组收集 <see cref="framesPerCase"/> 帧的
///    <c>XPBD.Cloth.Frame</c> 与 <c>XPBD.Cloth.Constraint</c> 两个 ProfilerMarker，
///    计算均值 / 中位数 / P99，最终写入 <see cref="outputCsvPath"/>。
/// 4. 建议在 IL2CPP + Development Build（Autoconnect Profiler 关闭）下运行，
///    Burst 必须开启，JobsDebugger / LeakDetection 关闭，否则数据不具备参考价值。
/// </summary>
[DefaultExecutionOrder(-10000)]
public class BenchmarkRunner : MonoBehaviour
{
    // ---------------- Inspector 配置 ----------------

    [Header("被测 Spawner（模板）")]
    [Tooltip("场景中的 ClothRuntimeSpawner 引用。Benchmark 会 Destroy 它并按配置重新 Instantiate。")]
    public ClothRuntimeSpawner targetSpawner;

    [Header("采样策略")]
    [Tooltip("每个配置的总帧数（建议 ≥ 600，即 10 秒@60fps）。")]
    [Min(60)] public int framesPerCase = 600;

    [Tooltip("每个配置的预热帧数，不计入统计，避免 Burst JIT / 缓存冷启动干扰。")]
    [Min(0)] public int warmupFrames = 120;

    [Tooltip("每个配置之间的间隔帧数（等 Dispose/GC 稳定）。")]
    [Min(0)] public int coolDownFrames = 30;

    [Header("实验 A：布料分辨率扫描")]
    public bool runResolutionSweep = true;
    [Tooltip("每档的 (segments, subdivision) 对。默认对应 32² / 64² / 128² / 256² 粒子规模。")]
    public Vector2Int[] resolutionCases = new Vector2Int[]
    {
        new Vector2Int(31, 31),
        new Vector2Int(63, 63),
        new Vector2Int(127, 127),
        new Vector2Int(255, 255),
    };

    [Header("实验 B：图着色 vs 串行 IJob")]
    public bool runColoringCompare = true;
    [Tooltip("对比时固定使用的 (segments, subdivision)。")]
    public Vector2Int coloringFixedResolution = new Vector2Int(127, 127);

    [Header("实验 C：多核扩展性")]
    public bool runThreadScaling = true;
    [Tooltip("线程数扫描（受 JobWorkerMaximumCount 上限限制）。")]
    public int[] threadCases = new[] { 1, 2, 4, 6, 8, 12, 16 };
    [Tooltip("多核实验固定的 (segments, subdivision)。建议中等规模，避免内存带宽瓶颈。")]
    public Vector2Int threadFixedResolution = new Vector2Int(127, 127);

    [Header("输出")]
    [Tooltip("相对 Application.persistentDataPath 的 CSV 文件名。")]
    public string outputCsvPath = "xpbd_benchmark.csv";

    // ---------------- 运行时状态 ----------------

    private readonly StringBuilder _csv = new StringBuilder();
    private Recorder _frameRecorder;
    private Recorder _constraintRecorder;

    /// <summary>
    /// 模板 GameObject：Inspector 拖进来的 Spawner 所挂载的 GameObject，始终保持 SetActive(false)，
    /// 只用于 Instantiate 副本，绝不销毁、绝不激活。
    /// </summary>
    private GameObject _templateGO;

    /// <summary>模板默认配置的快照，避免在 _templateGO 被修改后丢失原始参数。</summary>
    private int _tplSegments, _tplSubdivision;
    private bool _tplUseColoring;

    /// <summary>当前正在跑 benchmark 的活动副本（会被反复 Destroy / Instantiate）。</summary>
    private ClothRuntimeSpawner _activeSpawner;

    IEnumerator Start()
    {
        if (targetSpawner == null)
        {
            Debug.LogError("[Bench] targetSpawner 未设置，无法启动。");
            yield break;
        }

        // 把 Inspector 拖进来的那个对象当作"模板"，立即冻住：
        //  - SetActive(false) 让它自己不会触发 Start/Spawn；
        //  - 保存配置快照，后续从快照复制而不是从模板对象读（即使模板被改也不影响）。
        _templateGO = targetSpawner.gameObject;
        _tplSegments = targetSpawner.segments;
        _tplSubdivision = targetSpawner.subdivision;
        _tplUseColoring = targetSpawner.useGraphColoring;
        _templateGO.SetActive(false);
        _activeSpawner = null;   // 还没有活动副本
        targetSpawner = null;    // 防止误用

        // === 关键：解除帧率锁定，让測量值不被 VSync/默认 targetFrameRate 削平 ===
        // 否则 fps 会锁在 60/120、frame_ms 倒推回去全是 16.6ms，完全反映不出系统真实负载。
        QualitySettings.vSyncCount = 0;
        Application.targetFrameRate = -1;   // 不限帧

        // 启用 Benchmark 同步模式：让 ClothSimulationSystem 在 Marker 结束前强制 Complete Job，
        // 保证 constraint_ms / frame_ms 測量的是真实执行时间，而不是仅调度时间。
        ClothSimulationSystem.BenchmarkSyncMode = true;

        // 表头
        _csv.AppendLine("experiment,case,segments,subdivision,numParticles,useColoring,threads," +
                        "frame_ms_mean,frame_ms_p50,frame_ms_p99," +
                        "constraint_ms_mean,constraint_ms_p50,constraint_ms_p99," +
                        "fps_mean");

        // ProfilerMarker 名称必须与 ClothSimulationSystem 中保持一致
        _frameRecorder = Recorder.Get("XPBD.Cloth.Frame");
        _constraintRecorder = Recorder.Get("XPBD.Cloth.Constraint");
        _frameRecorder.enabled = true;
        _constraintRecorder.enabled = true;

        // 预热：让 Burst 编译并填充 CPU 缓存
        yield return SpawnAndWait(_tplSegments, _tplSubdivision,
                                  _tplUseColoring, JobsUtility.JobWorkerCount,
                                  warmupFrames + 30, collect: false, experiment: "warmup");

        // === 实验 A：分辨率扫描 ===
        // 为了与实验 B（图着色对比）交叉验证，每档分辨率都跑两遍：
        //   - useColoring=false：作为基线，反映单线程 IJob 的朴素实现
        //   - useColoring=true ：本文最终方案，反映启用并行图着色后的性能
        // 输出 CSV 里通过 useColoring 列区分两组数据，论文表 5-1 直接按该列分两栏渲染。
        if (runResolutionSweep)
        {
            foreach (var useColor in new[] { false, true })
            {
                foreach (var r in resolutionCases)
                {
                    yield return SpawnAndWait(r.x, r.y, useColor,
                                              JobsUtility.JobWorkerCount,
                                              framesPerCase, collect: true,
                                              experiment: "resolution");
                    yield return CoolDown();
                }
            }
        }

        // === 实验 B ===
        if (runColoringCompare)
        {
            foreach (var useColor in new[] { false, true })
            {
                yield return SpawnAndWait(coloringFixedResolution.x, coloringFixedResolution.y,
                                          useColor, JobsUtility.JobWorkerCount,
                                          framesPerCase, collect: true,
                                          experiment: "coloring");
                yield return CoolDown();
            }
        }

        // === 实验 C ===
        if (runThreadScaling)
        {
            int maxWorkers = JobsUtility.JobWorkerMaximumCount;
            int originalWorkers = JobsUtility.JobWorkerCount;
            foreach (var t in threadCases)
            {
                int clamped = Mathf.Clamp(t, 1, maxWorkers);
                JobsUtility.JobWorkerCount = clamped;
                yield return SpawnAndWait(threadFixedResolution.x, threadFixedResolution.y,
                                          true, clamped,
                                          framesPerCase, collect: true,
                                          experiment: "threads");
                yield return CoolDown();
            }
            JobsUtility.JobWorkerCount = originalWorkers;
        }

        // 写文件
        var path = Path.Combine(Application.persistentDataPath, outputCsvPath);
        File.WriteAllText(path, _csv.ToString());
        Debug.Log($"[Bench] CSV 已写入: {path}");
    }

    /// <summary>按给定配置销毁旧实例、创建新实例，并采集 frameCount 帧。</summary>
    IEnumerator SpawnAndWait(int segments, int subdivision, bool useColoring, int threadsInUse,
                              int frameCount, bool collect, string experiment)
    {
        // 销毁上一个活动副本（它自己 OnDestroy 会清理 Entity）。
        // 注意：只销毁 _activeSpawner，绝不碰 _templateGO。
        if (_activeSpawner != null)
        {
            Destroy(_activeSpawner.gameObject);
            _activeSpawner = null;
            // 等两帧，让 OnDestroy + ECS 清理生效
            yield return null;
            yield return null;
        }

        // 从模板克隆一个新的 Spawner GameObject。
        // 模板一直处于 SetActive(false)，Instantiate 出来的副本默认也是 inactive，
        // 这样可以先改完 segments/subdivision 等参数再激活，避免 Start 用默认参数跑一遍。
        if (_templateGO == null)
        {
            Debug.LogError("[Bench] 模板 GameObject 已丢失，无法继续。");
            yield break;
        }
        var go = Instantiate(_templateGO);
        go.name = $"ClothBench_{segments}x{subdivision}_color{useColoring}";
        var sp = go.GetComponent<ClothRuntimeSpawner>();
        sp.segments = Mathf.Max(1, segments);
        sp.subdivision = Mathf.Max(1, subdivision);
        sp.useGraphColoring = useColoring;
        sp.enabled = true;
        go.SetActive(true);  // 这时才触发 Spawner.Start → SpawnCloth
        _activeSpawner = sp;

        // 等待 Spawn 完成（Start 会在下一帧调用）
        yield return null;
        yield return null;
        yield return null;

        // 预热若干帧（仅首次 / 指定实验会 collect=true）
        int warm = Mathf.Max(warmupFrames, 30);
        for (int i = 0; i < warm; i++) yield return null;

        if (!collect) yield break;

        // 采样
        var frameSamples = new List<double>(frameCount);
        var constraintSamples = new List<double>(frameCount);
        // CPU 端到端帧耗时（unscaledDeltaTime，单位秒），最后取中位数再换算 fps，比瞬时 fps 累加稳定。
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

        int numParticles = (segments + 1) * (subdivision + 1);
        var (fm, fp50, fp99) = Stats(frameSamples);
        var (cm, cp50, cp99) = Stats(constraintSamples);
        // FPS = 1 / dt_median，避免极值帧把均值拉偏
        var (_, dtMedian, _) = Stats(frameDtSamples);
        float fpsMean = dtMedian > 0 ? (float)(1.0 / dtMedian) : 0f;

        _csv.AppendLine($"{experiment},{segments}x{subdivision},{segments},{subdivision},{numParticles}," +
                        $"{useColoring},{threadsInUse}," +
                        $"{fm:F3},{fp50:F3},{fp99:F3},{cm:F3},{cp50:F3},{cp99:F3},{fpsMean:F1}");

        Debug.Log($"[Bench] {experiment} seg={segments} sub={subdivision} color={useColoring} " +
                  $"threads={threadsInUse} | frame {fm:F2}/{fp50:F2}/{fp99:F2} ms | " +
                  $"cons {cm:F2}/{cp50:F2}/{cp99:F2} ms | fps={fpsMean:F1}");
    }

    IEnumerator CoolDown()
    {
        for (int i = 0; i < coolDownFrames; i++) yield return null;
    }

    /// <summary>返回 (mean, p50, p99)。小样本时 p99 退化为 max。</summary>
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
