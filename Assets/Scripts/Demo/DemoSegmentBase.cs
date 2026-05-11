using UnityEngine;

namespace Ropest.Demo
{
    /// <summary>
    /// Demo 录屏专用：单个演示段抽象基类。
    /// 不依赖 MonoBehaviour——由 DemoDirector 用 new() 创建并轮播。
    /// 一个段必须给出：标题（HUD 字幕）、持续时间、相机预设、Enter/Update/Exit 三步生命周期。
    ///
    /// 设计原则：
    ///   - 段内创建的所有 GameObject 都挂在 DemoDirector.SegmentRoot 下；Exit 时统一 Destroy 父节点即可清理干净，
    ///     避免遗漏导致下一段画面叠加。
    ///   - 段内不引入对论文实现的任何修改；纯组合调用 ClothRuntimeSpawner / SoftBodyRuntimeSpawner / AnalyticalColliderSource。
    /// </summary>
    public abstract class DemoSegmentBase
    {
        /// <summary>HUD 显示的中文段标题。</summary>
        public abstract string Title { get; }

        /// <summary>HUD 显示的中文副标题（描述当前展示的现象/特性），可为空。</summary>
        public virtual string Subtitle => string.Empty;

        /// <summary>本段持续多少秒（实时秒，不是物理时间）。这是"正常渲染"阶段的时长，不含可视化阶段。</summary>
        public virtual float Duration => 10f;

        /// <summary>正常渲染结束后，进入"粒子可视化"阶段的时长（秒）。
        /// 默认等于 Duration（可视化阶段是整段重播，时长保持一致最自然）。
        /// 不需要可视化的段（如 Seg07 大规模性能段）覆写为 0。
        /// </summary>
        public virtual float VisualizationDuration => Duration;

        /// <summary>段总时长 = Duration + VisualizationDuration。Director 用它来排相机和切段。</summary>
        public float TotalDuration => Duration + VisualizationDuration;

        /// <summary>本段启用的相机预设。镜头插值在整个 TotalDuration 上进行（SmoothStep）。</summary>
        public abstract DemoCameraPreset Camera { get; }

        protected DemoDirector Director { get; private set; }
        protected Transform Root { get; private set; }
        protected float ElapsedTime { get; private set; }

        public void Bind(DemoDirector director, Transform root)
        {
            Director = director;
            Root = root;
            ElapsedTime = 0f;
        }

        public void TickEnter() { OnEnter(); }
        public void TickUpdate(float dt)
        {
            ElapsedTime += dt;
            OnUpdate(dt);
        }
        public void TickExit() { OnExit(); }

        protected abstract void OnEnter();
        protected virtual void OnUpdate(float dt) { }
        protected virtual void OnExit() { }
    }
}
