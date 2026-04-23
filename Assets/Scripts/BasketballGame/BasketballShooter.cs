using Unity.Entities;
using Unity.Mathematics;
using UnityEngine;

namespace BasketballGame
{
    /// <summary>
    /// 投篮输入控制 + 蓄力条 + 出手
    /// 
    /// 交互流程：
    /// 1. 按下鼠标左键 -> 进入蓄力状态，蓄力值 0→1→0 来回震荡（速度条式）
    /// 2. 松开鼠标左键 -> 按瞬时蓄力值投出一颗球
    /// 3. 出手方向 = 相机 forward 绕 cameraRight 向上抬起 pitchBoostDegree 度
    /// </summary>
    public class BasketballShooter : MonoBehaviour
    {
        [Header("引用")]
        [Tooltip("用于决定发射方向/起点的相机")]
        public Camera shootCamera;
        [Tooltip("球Prefab：必须挂有 SoftBodyRuntimeSpawner + BasketballBall")]
        public GameObject ballPrefab;

        [Header("按键")]
        public KeyCode chargeButton = KeyCode.Mouse0;

        [Header("蓄力")]
        [Range(0.2f, 3f)] public float cycleDuration = 1.2f;

        // -----------------------------
        // 运行时状态
        // -----------------------------
        private bool _charging;
        private float _chargeTimer;
        /// <summary>当前蓄力值（0~1），为蓄力条的显示值</summary>
        public float CurrentPower { get; private set; }
        public bool IsCharging => _charging;

        void Awake()
        {
            if (shootCamera == null)
            {
                shootCamera = Camera.main;
            }
        }

        void Update()
        {
            // 只在光标被锁定时（游戏模式下）响应投篮
            var cam = shootCamera != null ? shootCamera.GetComponent<BasketballCameraController>() : null;
            if (cam != null && !cam.IsCursorLocked) return;

            if (Input.GetKeyDown(chargeButton))
            {
                _charging = true;
                _chargeTimer = 0f;
                CurrentPower = 0f;
            }

            if (_charging)
            {
                _chargeTimer += Time.deltaTime;
                // PingPong 在 [0, 1] 之间往返，周期= 2*cycleDuration
                CurrentPower = Mathf.PingPong(_chargeTimer / cycleDuration, 1f);

                if (Input.GetKeyUp(chargeButton))
                {
                    Shoot(CurrentPower);
                    _charging = false;
                    CurrentPower = 0f;
                    _chargeTimer = 0f;
                }
            }
        }

        /// <summary>按蓄力比例投出一个球</summary>
        public void Shoot(float power01)
        {
            if (ballPrefab == null)
            {
                Debug.LogError("[BasketballShooter] ballPrefab 未设置");
                return;
            }

            var gm = BasketballGameManager.Instance;
            if (gm == null) return;

            // 生成位置（相机局部偏移）
            Vector3 spawnPos = shootCamera.transform.TransformPoint(gm.ballSpawnLocalOffset);

            // 生成方向：相机 forward 绕 cameraRight 上抬 pitchBoostDegree 度
            Vector3 fwd = shootCamera.transform.forward;
            Vector3 right = shootCamera.transform.right;
            Quaternion pitchUp = Quaternion.AngleAxis(-gm.pitchBoostDegree, right); // 绕右轴负向=抬头
            Vector3 dir = (pitchUp * fwd).normalized;

            // 速度
            float speed = Mathf.Lerp(gm.minShootSpeed, gm.maxShootSpeed, Mathf.Clamp01(power01));
            Vector3 initialVelocity = dir * speed;

            // 实例化（prefab 自带 SoftBodyRuntimeSpawner，但模板是 inactive，需要激活后才能触发 Start）
            var ballGo = Instantiate(ballPrefab, spawnPos, Quaternion.identity);
            var ball = ballGo.GetComponent<BasketballBall>();
            if (ball == null)
            {
                ball = ballGo.AddComponent<BasketballBall>();
            }
            ball.SetInitialVelocity(initialVelocity);
            ballGo.SetActive(true);

            gm.RegisterBall(ball);
        }
    }
}
