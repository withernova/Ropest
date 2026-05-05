using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Text;
using UnityEngine;

/// <summary>
/// OOP 基线布料 Benchmark（对应论文表 5-5：DOTS vs OOP 对比）。
///
/// 驱动对象为 <see cref="global::Cloth"/>（<c>Assets/Scripts/XPBDRelated/Cloth/Cloth.cs</c>），
/// 即非 DOTS 的单线程 MonoBehaviour 版 XPBD 实现，算法逻辑与 DOTS 版本一致。
///
/// 用法：
/// 1. 在场景内放一个挂了 <c>Cloth</c> 组件的 GameObject（必须含 MeshFilter/MeshRenderer 依赖），
///    把它拖进 <see cref="targetCloth"/>；
/// 2. 配置 <see cref="resolutionCases"/> 为 (segments, subdivision) 列表；
///    默认跟 DOTS Benchmark 保持一致：32² / 64² / 128² / 256²；
/// 3. 运行场景，脚本会自动遍历所有配置，读取 Cloth.LastSolveMs / LastFrameMs，
///    最终写入 CSV；
/// 4. 采样协议（预热 / 窗口 / 冷却）与 <see cref="BenchmarkRunner"/> 完全一致。
///
/// 注意：
/// - OOP 基线在 FixedUpdate 内运行，默认 FixedDeltaTime = 0.02s (50Hz)。
///   为与 DOTS 版对齐，脚本会把 Time.fixedDeltaTime 设为 1/60 (16.67ms)，
///   保证子步频率一致；实际单帧耗时直接由 Stopwatch 量（不受固定时长影响）。
/// - 256² 规模的单线程 OOP 版会非常慢（估计 >200 ms/帧），脚本会自动把 framesPerCase
///   减半到 300，避免一次跑 10 分钟。
/// </summary>
[DefaultExecutionOrder(-10000)]
public class BaselineBenchmarkRunner : MonoBehaviour
{
    // ---------------- Inspector 配置 ----------------

    [Header("被测 Cloth（模板）")]
    [Tooltip("场景中挂了 Cloth 组件的 GameObject。会被 Destroy 并按配置重新 Instantiate。")]
    public Cloth targetCloth;

    [Header("采样策略")]
    [Min(60)] public int framesPerCase = 600;
    [Min(0)] public int warmupFrames = 120;
    [Min(0)] public int coolDownFrames = 30;

    [Header("分辨率扫描")]
    [Tooltip("(segments, subdivision) 对。与 DOTS Benchmark 对齐：31/63/127/255（对应 32²~256² 粒子）。")]
    public Vector2Int[] resolutionCases = new Vector2Int[]
    {
        new Vector2Int(31, 31),
        new Vector2Int(63, 63),
        new Vector2Int(127, 127),
        new Vector2Int(255, 255),
    };

    [Header("输出")]
    public string outputCsvPath = "xpbd_baseline_benchmark.csv";

    // ---------------- 运行时状态 ----------------

    private readonly StringBuilder _csv = new StringBuilder();
    private GameObject _templateGO;
    private int _tplSegments, _tplSubdivision;
    private Cloth _activeCloth;

    IEnumerator Start()
    {
        if (targetCloth == null)
        {
            Debug.LogError("[BaseBench] targetCloth 未设置，无法启动。");
            yield break;
        }

        _templateGO = targetCloth.gameObject;
        _tplSegments = targetCloth.segments;
        _tplSubdivision = targetCloth.subdivision;
        _templateGO.SetActive(false);
        _activeCloth = null;
        targetCloth = null;

        QualitySettings.vSyncCount = 0;
        Application.targetFrameRate = -1;

        // 与 DOTS 版对齐：FixedStep 60Hz
        Time.fixedDeltaTime = 1f / 60f;

        _csv.AppendLine("experiment,case,segments,subdivision,numParticles," +
                        "frame_ms_mean,frame_ms_p50,frame_ms_p99," +
                        "solve_ms_mean,solve_ms_p50,solve_ms_p99," +
                        "fps_mean");

        foreach (var r in resolutionCases)
        {
            // 256² 规模缩短窗口避免卡住（单线程 XPBD 预计每帧 >200 ms）
            int frames = (r.x >= 255) ? Mathf.Min(framesPerCase, 300) : framesPerCase;
            yield return SpawnAndWait(r.x, r.y, frames, collect: true, experiment: "baseline_oop");
            yield return CoolDown();
        }

        var path = Path.Combine(Application.persistentDataPath, outputCsvPath);
        File.WriteAllText(path, _csv.ToString());
        Debug.Log($"[BaseBench] CSV 已写入: {path}");
    }

    IEnumerator SpawnAndWait(int segments, int subdivision, int frameCount, bool collect, string experiment)
    {
        if (_activeCloth != null)
        {
            Destroy(_activeCloth.gameObject);
            _activeCloth = null;
            yield return null;
            yield return null;
        }

        if (_templateGO == null)
        {
            Debug.LogError("[BaseBench] 模板 GameObject 已丢失。");
            yield break;
        }

        var go = Instantiate(_templateGO);
        go.name = $"ClothBaseline_{segments}x{subdivision}";
        var c = go.GetComponent<Cloth>();
        c.segments = Mathf.Max(1, segments);
        c.subdivision = Mathf.Max(1, subdivision);
        c.simulate = true;
        c.enabled = true;
        go.SetActive(true);  // 激活后触发 Awake → CreateCloth → InitSolver
        _activeCloth = c;

        // 等 solver 初始化就绪（InitSolver 是协程，可能跨多帧）
        for (int i = 0; i < 10; i++) yield return new WaitForFixedUpdate();

        // 预热
        int warm = Mathf.Max(warmupFrames, 30);
        for (int i = 0; i < warm; i++) yield return new WaitForFixedUpdate();

        if (!collect) yield break;

        var frameSamples = new List<double>(frameCount);
        var solveSamples = new List<double>(frameCount);
        var dtSamples = new List<double>(frameCount);
        double lastStamp = Time.realtimeSinceStartupAsDouble;

        for (int i = 0; i < frameCount; i++)
        {
            yield return new WaitForFixedUpdate();

            // 从 Cloth 组件上读取插桩值
            if (c != null)
            {
                if (c.LastFrameMs > 0) frameSamples.Add(c.LastFrameMs);
                if (c.LastSolveMs > 0) solveSamples.Add(c.LastSolveMs);
            }

            // 用 realtime 量相邻 FixedUpdate 的真实间隔 → 倒推"端到端物理更新频率"
            double now = Time.realtimeSinceStartupAsDouble;
            double dt = now - lastStamp;
            lastStamp = now;
            if (dt > 0) dtSamples.Add(dt);
        }

        int numParticles = (segments + 1) * (subdivision + 1);
        var (fm, fp50, fp99) = Stats(frameSamples);
        var (sm, sp50, sp99) = Stats(solveSamples);
        var (_, dtMedian, _) = Stats(dtSamples);
        float fpsMean = dtMedian > 0 ? (float)(1.0 / dtMedian) : 0f;

        _csv.AppendLine($"{experiment},{segments}x{subdivision},{segments},{subdivision},{numParticles}," +
                        $"{fm:F3},{fp50:F3},{fp99:F3},{sm:F3},{sp50:F3},{sp99:F3},{fpsMean:F1}");

        Debug.Log($"[BaseBench] {experiment} seg={segments} sub={subdivision} | " +
                  $"frame {fm:F2}/{fp50:F2}/{fp99:F2} ms | solve {sm:F2}/{sp50:F2}/{sp99:F2} ms | fps={fpsMean:F1}");
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
