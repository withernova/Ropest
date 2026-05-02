using UnityEngine;
using UnityEngine.UI;

namespace BasketballGame
{
    public class BasketballHUD : MonoBehaviour
    {
        [Header("引用")]
        public BasketballShooter shooter;

        // UI 元素（运行时创建）
        private Canvas _canvas;
        private Slider _powerBar;
        private Image _powerFill;
        private Text _scoreText;
        private Image _crosshair;
        private Slider _stayBar;
        private Text _stateText;
        private Text _goalPopupText;
        private int _lastScore = 0;    // 用于检测分数增量触发动画
        private float _goalPopupTimer = 0f;  // 剩余显示时间（0=不显示）
        private float _scorePulseTimer = 0f; // 分数文字放大脉冲剩余时间

        void Start()
        {
            BuildUI();
            var gm = BasketballGameManager.Instance;
            if (gm != null)
            {
                gm.OnScoreChanged += OnScoreChanged;
                OnScoreChanged(gm.Score);
            }
        }

        void OnDestroy()
        {
            var gm = BasketballGameManager.Instance;
            if (gm != null) gm.OnScoreChanged -= OnScoreChanged;
        }

        void Update()
        {
            // 蓄力条
            if (_powerBar != null && shooter != null)
            {
                _powerBar.value = shooter.CurrentPower;
                // 颜色：低=绿，中=黄，高=红
                if (_powerFill != null)
                {
                    _powerFill.color = Color.Lerp(Color.green, Color.red, shooter.CurrentPower);
                }
                _powerBar.gameObject.SetActive(shooter.IsCharging);
            }

            // 取最"危险"的那张布（停留时间最高的）显示进度
            var gm = BasketballGameManager.Instance;
            if (gm != null && _stayBar != null)
            {
                float maxStay = 0f;
                foreach (var c in gm.GetActiveCloths())
                {
                    if (c == null) continue;
                    if (c.StayProgress01 > maxStay) maxStay = c.StayProgress01;
                }
                _stayBar.value = maxStay;
                _stayBar.gameObject.SetActive(maxStay > 0.01f);
            }

            // "GOAL! +1" 弹窗淡出动画
            if (_goalPopupText != null && _goalPopupTimer > 0f)
            {
                _goalPopupTimer -= Time.deltaTime;
                float t = Mathf.Clamp01(_goalPopupTimer / 1.5f); // 1→0 随时间衰减
                // 淡出 + 轻微上浮 + 头段放大
                float scale = Mathf.Lerp(0.8f, 1.4f, Mathf.Clamp01((1.5f - _goalPopupTimer) / 0.25f));
                if (_goalPopupTimer < 1.0f) scale = Mathf.Lerp(1.0f, 1.4f, t); // 后段回落到 1.0
                _goalPopupText.rectTransform.localScale = Vector3.one * scale;
                var c0 = _goalPopupText.color;
                c0.a = Mathf.Clamp01(t * 2f); // 最后 0.75s 淡出
                _goalPopupText.color = c0;
                var rt = _goalPopupText.rectTransform;
                rt.anchoredPosition = new Vector2(0, 40f + (1f - t) * 30f); // 轻微上浮
                if (_goalPopupTimer <= 0f)
                {
                    _goalPopupText.gameObject.SetActive(false);
                }
            }

            // 分数文字放大脉冲
            if (_scoreText != null && _scorePulseTimer > 0f)
            {
                _scorePulseTimer -= Time.deltaTime;
                float t = Mathf.Clamp01(_scorePulseTimer / 0.6f);
                // 前 0.15s 快速放大，后 0.45s 平滑回落
                float scale = 1f + Mathf.Sin(Mathf.Clamp01(t) * Mathf.PI) * 0.35f;
                _scoreText.rectTransform.localScale = Vector3.one * scale;
                if (_scorePulseTimer <= 0f)
                {
                    _scoreText.rectTransform.localScale = Vector3.one;
                }
            }
        }

        private void OnScoreChanged(int s)
        {
            if (_scoreText != null) _scoreText.text = $"Score: {s}";

            // 分数增加时触发 "GOAL! +delta" 弹窗 + 分数文字放大脉冲
            int delta = s - _lastScore;
            if (delta > 0)
            {
                if (_goalPopupText != null)
                {
                    _goalPopupText.text = delta > 1 ? $"GOAL! +{delta}" : "GOAL! +1";
                    _goalPopupText.gameObject.SetActive(true);
                }
                _goalPopupTimer = 1.5f;   // 提示显示 1.5 秒
                _scorePulseTimer = 0.6f;  // 分数文字脉冲 0.6 秒
            }
            _lastScore = s;
        }

        private void BuildUI()
        {
            // Canvas
            var canvasGo = new GameObject("BasketballHUDCanvas");
            canvasGo.transform.SetParent(transform, false);
            _canvas = canvasGo.AddComponent<Canvas>();
            _canvas.renderMode = RenderMode.ScreenSpaceOverlay;
            canvasGo.AddComponent<CanvasScaler>().uiScaleMode = CanvasScaler.ScaleMode.ScaleWithScreenSize;
            canvasGo.GetComponent<CanvasScaler>().referenceResolution = new Vector2(1920, 1080);
            canvasGo.AddComponent<GraphicRaycaster>();

            // 分数
            _scoreText = CreateText(canvasGo.transform, "Score: 0",
                new Vector2(1, 1), new Vector2(1, 1), new Vector2(-180, -50),
                52, TextAnchor.UpperRight);
            _scoreText.color = new Color(1f, 0.92f, 0.35f); // 金黄
            _scoreText.fontStyle = FontStyle.Bold;
            _scoreText.rectTransform.sizeDelta = new Vector2(520, 100);
            _scoreText.rectTransform.pivot = new Vector2(1, 1);

            // 准心
            var crossGo = new GameObject("Crosshair");
            crossGo.transform.SetParent(canvasGo.transform, false);
            _crosshair = crossGo.AddComponent<Image>();
            _crosshair.color = new Color(1f, 1f, 1f, 0.85f);
            var crt = _crosshair.rectTransform;
            crt.anchorMin = crt.anchorMax = new Vector2(0.5f, 0.5f);
            crt.sizeDelta = new Vector2(8, 8);
            crt.anchoredPosition = Vector2.zero;

            // 蓄力条（底部中央，水平）
            _powerBar = CreateSlider(canvasGo.transform, "PowerBar",
                new Vector2(0.5f, 0f), new Vector2(0.5f, 0f), new Vector2(0, 80),
                new Vector2(520, 28),
                out _powerFill);
            _powerBar.gameObject.SetActive(false);

            // 停留进度条（屏幕上方中央）
            _stayBar = CreateSlider(canvasGo.transform, "StayBar",
                new Vector2(0.5f, 1f), new Vector2(0.5f, 1f), new Vector2(0, -80),
                new Vector2(420, 18),
                out var stayFill);
            if (stayFill != null) stayFill.color = new Color(0.2f, 0.9f, 0.5f, 0.9f);
            _stayBar.gameObject.SetActive(false);

            // 状态提示
            _stateText = CreateText(canvasGo.transform, "按住鼠标左键蓄力，松开投出",
                new Vector2(0.5f, 0f), new Vector2(0.5f, 0f), new Vector2(0, 130),
                22, TextAnchor.LowerCenter);

            // 进球弹窗（屏幕中央偏上，默认隐藏）
            _goalPopupText = CreateText(canvasGo.transform, "GOAL! +1",
                new Vector2(0.5f, 0.5f), new Vector2(0.5f, 0.5f), new Vector2(0, 40),
                88, TextAnchor.MiddleCenter);
            _goalPopupText.color = new Color(1f, 0.85f, 0.2f);
            _goalPopupText.fontStyle = FontStyle.Bold;
            _goalPopupText.rectTransform.sizeDelta = new Vector2(800, 140);
            _goalPopupText.gameObject.SetActive(false);
        }

        private Text CreateText(Transform parent, string content,
            Vector2 anchorMin, Vector2 anchorMax, Vector2 anchoredPos,
            int fontSize, TextAnchor alignment)
        {
            var go = new GameObject("Text");
            go.transform.SetParent(parent, false);
            var t = go.AddComponent<Text>();
            t.text = content;
            t.font = Resources.GetBuiltinResource<Font>("LegacyRuntime.ttf");
            t.fontSize = fontSize;
            t.alignment = alignment;
            t.color = Color.white;
            t.raycastTarget = false;
            var rt = t.rectTransform;
            rt.anchorMin = anchorMin;
            rt.anchorMax = anchorMax;
            rt.anchoredPosition = anchoredPos;
            rt.sizeDelta = new Vector2(360, 60);
            // 阴影
            var sh = go.AddComponent<Shadow>();
            sh.effectColor = new Color(0, 0, 0, 0.8f);
            sh.effectDistance = new Vector2(2, -2);
            return t;
        }

        private Slider CreateSlider(Transform parent, string name,
            Vector2 anchorMin, Vector2 anchorMax, Vector2 anchoredPos, Vector2 size,
            out Image fillImage)
        {
            var go = new GameObject(name);
            go.transform.SetParent(parent, false);
            var rt = go.AddComponent<RectTransform>();
            rt.anchorMin = anchorMin;
            rt.anchorMax = anchorMax;
            rt.sizeDelta = size;
            rt.anchoredPosition = anchoredPos;

            var slider = go.AddComponent<Slider>();
            slider.minValue = 0f;
            slider.maxValue = 1f;
            slider.direction = Slider.Direction.LeftToRight;
            slider.transition = Selectable.Transition.None;
            slider.interactable = false;

            // Background
            var bgGo = new GameObject("Background");
            bgGo.transform.SetParent(go.transform, false);
            var bg = bgGo.AddComponent<Image>();
            bg.color = new Color(0, 0, 0, 0.55f);
            var bgRt = bg.rectTransform;
            bgRt.anchorMin = Vector2.zero; bgRt.anchorMax = Vector2.one;
            bgRt.offsetMin = Vector2.zero; bgRt.offsetMax = Vector2.zero;

            // Fill Area
            var fillArea = new GameObject("Fill Area");
            fillArea.transform.SetParent(go.transform, false);
            var fillAreaRt = fillArea.AddComponent<RectTransform>();
            fillAreaRt.anchorMin = new Vector2(0, 0);
            fillAreaRt.anchorMax = new Vector2(1, 1);
            fillAreaRt.offsetMin = new Vector2(4, 4);
            fillAreaRt.offsetMax = new Vector2(-4, -4);

            var fillGo = new GameObject("Fill");
            fillGo.transform.SetParent(fillArea.transform, false);
            fillImage = fillGo.AddComponent<Image>();
            fillImage.color = Color.green;
            var fillRt = fillImage.rectTransform;
            fillRt.anchorMin = Vector2.zero; fillRt.anchorMax = new Vector2(1, 1);
            fillRt.offsetMin = Vector2.zero; fillRt.offsetMax = Vector2.zero;

            slider.fillRect = fillRt;
            slider.targetGraphic = fillImage;
            slider.value = 0f;
            return slider;
        }
    }
}
