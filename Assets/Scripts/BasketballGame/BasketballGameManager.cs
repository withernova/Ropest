using System.Collections.Generic;
using UnityEngine;

namespace BasketballGame
{
    public class BasketballGameManager : MonoBehaviour
    {
        public static BasketballGameManager Instance { get; private set; }

        [Header("场地参数")]
        [Tooltip("玩家站位（相机以此为中心）")]
        public Transform playerStand;
        [Tooltip("布料升起区域中心（通常在玩家正前方）")]
        public Transform targetAreaCenter;
        [Tooltip("布料升起区域半径（随机位置的散布范围，XZ平面）")]
        public float targetAreaRadius = 6f;
        [Tooltip("布料升起的高度下限（相对 targetAreaCenter.y）")]
        public float clothRiseHeightMin = 3f;
        [Tooltip("布料升起的高度上限（相对 targetAreaCenter.y）。每张布料在 [Min, Max] 间随机一个高度")]
        public float clothRiseHeightMax = 7f;

        /// <summary>
        /// 兼容旧接口：返回一个随机高度，位于 [clothRiseHeightMin, clothRiseHeightMax] 之间。
        /// 每次读取都会返回一个新的随机值，用于让每张布料最终停在不同高度。
        /// </summary>
        public float clothRiseHeight
        {
            get
            {
                float lo = Mathf.Min(clothRiseHeightMin, clothRiseHeightMax);
                float hi = Mathf.Max(clothRiseHeightMin, clothRiseHeightMax);
                return Random.Range(lo, hi);
            }
        }

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

        [Header("进球判定")]
        [Tooltip("得分需要在布料内停留的时间（秒）")]
        public float requiredStayDuration = 3.5f;
        [Tooltip("判定球进入布料区域的球心半径（以布料中心为准）")]
        public float scoreTriggerRadius = 2.0f;

        [Header("数量上限")]
        public int maxActiveBalls = 2;
        public int maxActiveCloths = 3;

        private readonly List<BasketballBall> _activeBalls = new();
        private readonly List<BasketballClothTarget> _activeCloths = new();

        public int Score { get; private set; }

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

        public void AddScore(int delta, BasketballClothTarget source)
        {
            Score += delta;
            OnScoreChanged?.Invoke(Score);
            OnTargetScored?.Invoke(source);
        }
    }
}
