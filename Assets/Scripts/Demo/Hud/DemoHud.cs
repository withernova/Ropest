using Unity.Collections;
using Unity.Entities;
using UnityEngine;

namespace Ropest.Demo
{
    /// <summary>
    /// 录屏用 HUD（IMGUI，零依赖）：
    ///   - 顶部居中：当前段中文标题
    ///   - 顶部居中（小字）：副标题（描述当前展示现象）
    ///   - 左下：粒子总数（布料 + 软体）、FPS、当前 Demo 名
    ///   - 右下：进度条（本段时间）
    ///   - 右上小角标：项目水印 "DOTS-XPBD Demo"
    ///
    /// 之所以用 IMGUI 而不是 UGUI Canvas：盲审视频不需要复杂排版，IMGUI 最少代码、最稳，
    /// 而且 OnGUI 不会被 UGUI 的 EventSystem 遮挡或拦截。
    /// </summary>
    [DisallowMultipleComponent]
    public class DemoHud : MonoBehaviour
    {
        public DemoDirector Director;

        [Header("水印")]
        public string Watermark = "DOTS-XPBD Demo";

        [Header("样式")]
        public int TitleFontSize = 36;
        public int SubtitleFontSize = 20;
        public int InfoFontSize = 18;
        public Color TitleColor = new Color(1f, 1f, 1f, 0.95f);
        public Color SubtitleColor = new Color(1f, 1f, 1f, 0.7f);
        public Color InfoColor = new Color(0.85f, 1f, 0.85f, 0.95f);

        // FPS 平滑
        float _fpsAccum;
        int _fpsFrames;
        float _fpsTimer;
        float _fpsSmoothed;

        // 缓存粒子统计（每 0.25s 刷一次，不用每帧查 ECS）
        float _particleQueryTimer;
        int _cachedClothParticles;
        int _cachedSoftBodyParticles;

        void Update()
        {
            // FPS 1Hz 更新
            _fpsAccum += Time.unscaledDeltaTime;
            _fpsFrames++;
            _fpsTimer += Time.unscaledDeltaTime;
            if (_fpsTimer >= 0.5f)
            {
                _fpsSmoothed = _fpsFrames / _fpsAccum;
                _fpsAccum = 0f;
                _fpsFrames = 0;
                _fpsTimer = 0f;
            }

            // 粒子统计 4Hz 更新
            _particleQueryTimer += Time.unscaledDeltaTime;
            if (_particleQueryTimer >= 0.25f)
            {
                _particleQueryTimer = 0f;
                RefreshParticleCounts();
            }
        }

        void RefreshParticleCounts()
        {
            var world = World.DefaultGameObjectInjectionWorld;
            if (world == null || !world.IsCreated)
            {
                _cachedClothParticles = 0;
                _cachedSoftBodyParticles = 0;
                return;
            }
            var em = world.EntityManager;

            _cachedClothParticles = SumBufferLengths<ClothTag>(em);
            _cachedSoftBodyParticles = SumBufferLengths<SoftBodyTag>(em);
        }

        static int SumBufferLengths<TTag>(EntityManager em) where TTag : unmanaged, IComponentData
        {
            var query = em.CreateEntityQuery(
                ComponentType.ReadOnly<TTag>(),
                ComponentType.ReadOnly<ParticlePosition>());
            var entities = query.ToEntityArray(Allocator.Temp);
            int total = 0;
            for (int i = 0; i < entities.Length; i++)
            {
                total += em.GetBuffer<ParticlePosition>(entities[i], true).Length;
            }
            entities.Dispose();
            query.Dispose();
            return total;
        }

        void OnGUI()
        {
            if (Director == null) return;

            float w = Screen.width;
            float h = Screen.height;

            // ====== 顶部居中：段标题 ======
            string title = Director.CurrentSegmentTitle;
            string subtitle = Director.CurrentSegmentSubtitle;

            var titleStyle = new GUIStyle(GUI.skin.label)
            {
                fontSize = TitleFontSize,
                alignment = TextAnchor.MiddleCenter,
                fontStyle = FontStyle.Bold,
                normal = { textColor = TitleColor }
            };
            var subtitleStyle = new GUIStyle(GUI.skin.label)
            {
                fontSize = SubtitleFontSize,
                alignment = TextAnchor.MiddleCenter,
                normal = { textColor = SubtitleColor }
            };

            DrawTextWithShadow(new Rect(0, 30, w, TitleFontSize + 8), title, titleStyle);
            if (!string.IsNullOrEmpty(subtitle))
            {
                DrawTextWithShadow(new Rect(0, 30 + TitleFontSize + 6, w, SubtitleFontSize + 8),
                    subtitle, subtitleStyle);
            }

            // ====== 左下：粒子数 / FPS / Demo 名 ======
            var infoStyle = new GUIStyle(GUI.skin.label)
            {
                fontSize = InfoFontSize,
                alignment = TextAnchor.LowerLeft,
                normal = { textColor = InfoColor }
            };
            int totalParticles = _cachedClothParticles + _cachedSoftBodyParticles;
            string info =
                $"{Watermark}\n" +
                $"Cloth:    {_cachedClothParticles,7} pts\n" +
                $"SoftBody: {_cachedSoftBodyParticles,7} pts\n" +
                $"TOTAL:    {totalParticles,7} pts\n" +
                $"FPS:      {_fpsSmoothed,7:F1}";
            DrawTextWithShadow(new Rect(20, h - 160, 400, 150), info, infoStyle);

            // ====== 右下：本段进度条 ======
            float pad = 20f;
            float barW = 320f;
            float barH = 6f;
            var barRect = new Rect(w - pad - barW, h - pad - barH - 18, barW, barH);
            var bgColor = new Color(1, 1, 1, 0.18f);
            var fgColor = new Color(0.4f, 1f, 0.55f, 0.9f);
            DrawRect(barRect, bgColor);
            float progress = Mathf.Clamp01(Director.CurrentSegmentProgress);
            DrawRect(new Rect(barRect.x, barRect.y, barRect.width * progress, barRect.height), fgColor);

            var segLabelStyle = new GUIStyle(GUI.skin.label)
            {
                fontSize = 14,
                alignment = TextAnchor.UpperRight,
                normal = { textColor = new Color(1, 1, 1, 0.7f) }
            };
            string segText = $"段 {Director.CurrentSegmentIndex + 1}/{Director.TotalSegments}    {Director.CurrentSegmentRemaining:F1}s";
            GUI.Label(new Rect(w - pad - barW, h - pad - barH - 38, barW, 18), segText, segLabelStyle);
        }

        static Texture2D _whiteTex;
        static Texture2D WhiteTex
        {
            get
            {
                if (_whiteTex == null)
                {
                    _whiteTex = new Texture2D(1, 1);
                    _whiteTex.SetPixel(0, 0, Color.white);
                    _whiteTex.Apply();
                }
                return _whiteTex;
            }
        }

        static void DrawRect(Rect r, Color c)
        {
            var prev = GUI.color;
            GUI.color = c;
            GUI.DrawTexture(r, WhiteTex);
            GUI.color = prev;
        }

        static void DrawTextWithShadow(Rect r, string text, GUIStyle style)
        {
            // 先画一层黑色阴影错开 1px，再画原色，简单但有效地提高字幕在任意背景下的可读性。
            var shadowStyle = new GUIStyle(style);
            shadowStyle.normal.textColor = new Color(0, 0, 0, 0.6f);
            var shadowRect = new Rect(r.x + 2, r.y + 2, r.width, r.height);
            GUI.Label(shadowRect, text, shadowStyle);
            GUI.Label(r, text, style);
        }
    }
}
