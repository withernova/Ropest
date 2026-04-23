using System.Collections.Generic;
using UnityEngine;

namespace BasketballGame
{
    /// <summary>
    /// 投篮小游戏总控
    /// 职责：
    /// 1. 管理场上球 / 布料 的数量上限（最多 2 球 + 3 布）
    /// 2. 维护计分
    /// 3. 提供全局配置（投掷参数、布料参数、场地参数）
    /// </summary>
    public class BasketballGameManager : MonoBehaviour
    {
        public static BasketballGameManager Instance { get; private set; }

        // -----------------------------
        // 场地与目标区域配置
        // -----------------------------
        [Header("场地参数")]
        [Tooltip("玩家站位（相机以此为中心）")]
        public Transform playerStand;
        [Tooltip("布料升起区域中心（通常在玩家正前方）")]
        public Transform targetAreaCenter;
        [Tooltip("布料升起区域半径（随机位置的散布范围，XZ平面）")]
        public float targetAreaRadius = 6f;
        [Tooltip("布料升起的高度（相对 targetAreaCenter.y）")]
        public float clothRiseHeight = 3.5f;

        // -----------------------------
        // 投掷参数
        // -----------------------------
        [Header("投掷参数")]
        [Tooltip("最小出手速度（蓄力=0时）")]
        public float minShootSpeed = 4f;
        [Tooltip("最大出手速度（蓄力=1时）")]
        public float maxShootSpeed = 18f;
        [Tooltip("投掷方向相对相机forward的上抬角度（度）")]
        public float pitchBoostDegree = 15f;
        [Tooltip("蓄力条来回震荡一次的时间（秒）")]
        public float powerBarCycleDuration = 1.2f;
        [Tooltip("球初始位置（相对相机的偏移，局部空间）")]
        public Vector3 ballSpawnLocalOffset = new Vector3(0.4f, -0.3f, 0.8f);

        // -----------------------------
        // 进球判定
        // -----------------------------
        [Header("进球判定")]
        [Tooltip("得分需要在布料内停留的时间（秒）")]
        public float requiredStayDuration = 3.5f;
        [Tooltip("判定球进入布料区域的球心半径（以布料中心为准）")]
        public float scoreTriggerRadius = 2.0f;

        // -----------------------------
        // 数量上限
        // -----------------------------
        [Header("数量上限")]
        public int maxActiveBalls = 2;
        public int maxActiveCloths = 3;

        // -----------------------------
        // 内部状态
        // -----------------------------
        private readonly List<BasketballBall> _activeBalls = new();
        private readonly List<BasketballClothTarget> _activeCloths = new();

        public int Score { get; private set; }

        // -----------------------------
        // 事件
        // -----------------------------
        public System.Action<int> OnScoreChanged;
        public System.Action<BasketballClothTarget> OnTargetScored;

        void Awake()
        {
            if (Instance != null && Instance != this)
            {
                Destroy(gameObject);
                return;
            }
            Instance = this;
        }

        void OnDestroy()
        {
            if (Instance == this) Instance = null;
        }

        // =========================================================
        // 球管理
        // =========================================================
        public void RegisterBall(BasketballBall ball)
        {
            _activeBalls.Add(ball);
            EnforceBallLimit();
        }

        public void UnregisterBall(BasketballBall ball)
        {
            _activeBalls.Remove(ball);
        }

        private void EnforceBallLimit()
        {
            while (_activeBalls.Count > maxActiveBalls)
            {
                var oldest = _activeBalls[0];
                _activeBalls.RemoveAt(0);
                if (oldest != null) Destroy(oldest.gameObject);
            }
        }

        // =========================================================
        // 布料管理
        // =========================================================
        public void RegisterCloth(BasketballClothTarget cloth)
        {
            _activeCloths.Add(cloth);
            EnforceClothLimit();
        }

        public void UnregisterCloth(BasketballClothTarget cloth)
        {
            _activeCloths.Remove(cloth);
        }

        public IReadOnlyList<BasketballClothTarget> GetActiveCloths() => _activeCloths;

        public int ActiveClothCount => _activeCloths.Count;

        private void EnforceClothLimit()
        {
            while (_activeCloths.Count > maxActiveCloths)
            {
                var oldest = _activeCloths[0];
                _activeCloths.RemoveAt(0);
                if (oldest != null) oldest.FallAndDestroy();
            }
        }

        // =========================================================
        // 计分
        // =========================================================
        public void AddScore(int delta, BasketballClothTarget source)
        {
            Score += delta;
            OnScoreChanged?.Invoke(Score);
            OnTargetScored?.Invoke(source);
        }
    }
}
