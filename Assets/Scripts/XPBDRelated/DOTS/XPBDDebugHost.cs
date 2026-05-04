using Unity.Collections;
using Unity.Entities;
using Unity.Mathematics;
using UnityEngine;

/// <summary>
/// XPBD 物理仿真的调试宿主 MonoBehaviour：
/// 1. 响应按键：P 暂停/恢复，N/右箭头 单步推进，]（可配）连续步进 K 次。
/// 2. 从 DOTS World 读取所有布料 / 软体实体的 ParticlePosition + ParticleVelocity，
///    在 Scene / Game 视图中用 Gizmos 画出每个逻辑顶点的速度方向和大小。
/// 3. 可选地在 OnGUI 中绘制极简状态面板，方便调试时一眼看到是否已暂停 / 当前迭代数。
/// 
/// 挂载方式：由 BasketballSceneBootstrap 在运行时 AddComponent。也可在场景里手动挂到任意 GameObject。
/// 注意：Gizmos 只在 Scene 视图绘制；若希望 Game 视图也能看到，可以启用 Gizmos 面板的 "Game" 开关。
/// </summary>
[DisallowMultipleComponent]
public class XPBDDebugHost : MonoBehaviour
{
    [Header("按键")]
    [Tooltip("按下切换仿真暂停状态")]
    public KeyCode pauseKey = KeyCode.P;

    [Tooltip("暂停状态下，按下推进一次子步（一个约束迭代）")]
    public KeyCode stepKey = KeyCode.N;

    [Tooltip("暂停状态下，按下连续推进 BatchStepCount 次子步")]
    public KeyCode batchStepKey = KeyCode.M;

    [Range(1, 50)]
    public int batchStepCount = 5;

    [Header("速度可视化")]
    public bool drawVelocityGizmos = true;

    [Tooltip("速度向量缩放：实际绘制长度 = velocity.magnitude * scale")]
    [Range(0.001f, 1f)]
    public float velocityGizmoScale = 0.08f;

    [Tooltip("小于该速度阈值的顶点不绘制（避免静止顶点堆满屏幕）")]
    [Range(0f, 2f)]
    public float minSpeedToDraw = 0.02f;

    public Color velocityArrowColor = new Color(0.1f, 1f, 0.4f, 0.9f);
    public Color velocityArrowHeadColor = new Color(1f, 0.9f, 0.2f, 1f);

    [Header("状态面板 (OnGUI)")]
    public bool showOverlay = true;

    [Tooltip("参考 NumSubSteps：用于显示\"本虚拟帧第几个迭代\"的分母。\n" +
             "单步模式下，每按一次 N 相当于正常帧的一个 SubStep，连按 NumSubSteps 次等于推进一帧。\n" +
             "此字段仅影响 UI 显示，请填与 SoftBodyRuntimeSpawner/ClothRuntimeSpawner 的 numSubSteps 一致。")]
    [Range(1, 50)]
    public int referenceNumSubSteps = 10;

    private void Update()
    {
        // 输入处理
        if (Input.GetKeyDown(pauseKey))
        {
            XPBDDebugController.TogglePause();
        }
        if (Input.GetKeyDown(stepKey))
        {
            if (!XPBDDebugController.IsPaused) XPBDDebugController.SetPaused(true);
            XPBDDebugController.RequestStep(1);
        }
        if (Input.GetKeyDown(batchStepKey))
        {
            if (!XPBDDebugController.IsPaused) XPBDDebugController.SetPaused(true);
            XPBDDebugController.RequestStep(batchStepCount);
        }

        // 同步 Gizmos 开关到静态控制器（为了让 OnDrawGizmos 使用）。
        XPBDDebugController.DrawVelocityGizmos = drawVelocityGizmos;
        XPBDDebugController.VelocityGizmoScale = velocityGizmoScale;
    }

    private void OnDrawGizmos()
    {
        if (!drawVelocityGizmos) return;
        if (!Application.isPlaying) return;

        var world = World.DefaultGameObjectInjectionWorld;
        if (world == null || !world.IsCreated) return;
        var em = world.EntityManager;

        // 先同步依赖，确保读取的 Buffer 不与未完成的 Job 冲突。
        // 布料 & 软体系统都把最后的 PostSolve 写入 ParticleVelocity/ParticlePosition，
        // Gizmos 读之前必须 Complete。
        em.CompleteAllTrackedJobs();

        DrawEntitiesOfTag<ClothTag>(em);
        DrawEntitiesOfTag<SoftBodyTag>(em);
    }

    private void DrawEntitiesOfTag<TTag>(EntityManager em) where TTag : unmanaged, IComponentData
    {
        var query = em.CreateEntityQuery(
            ComponentType.ReadOnly<TTag>(),
            ComponentType.ReadOnly<ParticlePosition>(),
            ComponentType.ReadOnly<ParticleVelocity>()
        );

        var entities = query.ToEntityArray(Allocator.Temp);
        for (int i = 0; i < entities.Length; i++)
        {
            var entity = entities[i];
            var pos = em.GetBuffer<ParticlePosition>(entity, true);
            var vel = em.GetBuffer<ParticleVelocity>(entity, true);
            int n = math.min(pos.Length, vel.Length);
            for (int p = 0; p < n; p++)
            {
                float3 v = vel[p].Value;
                float speed = math.length(v);
                if (speed < minSpeedToDraw) continue;

                float3 start = pos[p].Value;
                float3 end = start + v * velocityGizmoScale;

                Gizmos.color = velocityArrowColor;
                Gizmos.DrawLine((Vector3)start, (Vector3)end);

                // 箭头头部（简单的小十字，避免三角面片开销）
                Gizmos.color = velocityArrowHeadColor;
                float headSize = math.min(0.05f, speed * velocityGizmoScale * 0.2f);
                if (headSize > 1e-4f)
                {
                    Gizmos.DrawSphere((Vector3)end, headSize);
                }
            }
        }
        entities.Dispose();
        query.Dispose();
    }

    private void OnGUI()
    {
        if (!showOverlay) return;
        var rect = new Rect(10, 10, 440, 150);
        GUI.Box(rect, GUIContent.none);
        GUILayout.BeginArea(rect);

        var richStyle = new GUIStyle(GUI.skin.label) { richText = true };
        bool paused = XPBDDebugController.IsPaused;

        GUILayout.Label(
            $"[XPBD Debug] Paused = <b>{paused}</b>    PendingSteps = {XPBDDebugController.PendingSubSteps}",
            richStyle);

        // 迭代序号显示：
        //  - 暂停状态下：显示\"SingleStepIndex\"（自本次暂停按 P 以来已推进的 SubStep 数），
        //    并换算到\"本虚拟帧的第几个迭代\" = (SingleStepIndex-1) % N + 1 / N。
        //  - 未暂停状态下：显示累计逻辑 Tick 数。
        int n = Mathf.Max(1, referenceNumSubSteps);
        int stepIdx = XPBDDebugController.SingleStepIndex; // 1-based，暂停会话内
        if (paused && stepIdx > 0)
        {
            int frameIdx = (stepIdx - 1) / n;           // 已走完多少个\"虚拟帧\"
            int iterInFrame = ((stepIdx - 1) % n) + 1;  // 本虚拟帧内第几个迭代 (1..N)
            GUILayout.Label(
                $"<b>Iteration:</b> <color=#ffee66>{iterInFrame}/{n}</color> " +
                $"(frame #{frameIdx + 1} since pause)    <b>TotalStepInPause:</b> {stepIdx}",
                richStyle);
        }
        else if (paused)
        {
            GUILayout.Label(
                $"<b>Iteration:</b> <color=#aaaaaa>-- / {n}</color>   (press [{stepKey}] to advance)",
                richStyle);
        }
        else
        {
            GUILayout.Label(
                $"<b>Running</b>   NumSubSteps per fixed tick ≈ {n}    TotalTicks = {XPBDDebugController.TotalPhysicsTicks}",
                richStyle);
        }

        GUILayout.Label($"Keys: [{pauseKey}] Pause/Resume    [{stepKey}] +1 step    [{batchStepKey}] +{batchStepCount} steps");
        GUILayout.Label($"Gizmos velocity scale = {velocityGizmoScale:F3}    minSpeed = {minSpeedToDraw:F2}");
        GUILayout.EndArea();
    }
}
