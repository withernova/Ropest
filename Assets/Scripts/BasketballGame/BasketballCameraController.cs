using UnityEngine;

namespace BasketballGame
{
    /// <summary>
    /// 第一人称相机控制（FPS锁鼠标模式）
    /// - 游戏启动即锁光标、隐藏光标
    /// - 鼠标 X 控制 Yaw（绕 Y 轴），鼠标 Y 控制 Pitch（绕右轴，限制 [-80°, 80°]）
    /// - 按 Esc 释放光标（方便开发调试）
    /// </summary>
    [RequireComponent(typeof(Camera))]
    public class BasketballCameraController : MonoBehaviour
    {
        [Header("灵敏度")]
        [Range(0.1f, 10f)] public float mouseSensitivity = 2.0f;

        [Header("仰俯角限制")]
        public float minPitch = -80f;
        public float maxPitch = 80f;

        [Header("Y偏移（相对站位）")]
        [Tooltip("相机相对 playerStand 的高度偏移（W/S 可在 minCameraHeight~maxCameraHeight 范围内调整）")]
        public float cameraHeight = 1.6f;
        [Tooltip("W/S 按键调整相机高度的速度（米/秒）")]
        public float heightAdjustSpeed = 2.0f;
        [Tooltip("相机高度下限")]
        public float minCameraHeight = 0.6f;
        [Tooltip("相机高度上限")]
        public float maxCameraHeight = 3.5f;

        private float _yaw;
        private float _pitch;
        private bool _cursorLocked = true;

        public float Yaw => _yaw;
        public float Pitch => _pitch;

        void Start()
        {
            LockCursor(true);
            var euler = transform.eulerAngles;
            _yaw = euler.y;
            _pitch = euler.x;
            if (_pitch > 180f) _pitch -= 360f;
        }

        void Update()
        {
            // Esc 释放光标（方便编辑器调试），再次点击窗口左键会重新锁定
            if (Input.GetKeyDown(KeyCode.Escape))
            {
                LockCursor(false);
            }
            if (!_cursorLocked && Input.GetMouseButtonDown(0))
            {
                LockCursor(true);
            }

            if (!_cursorLocked) return;

            float mx = Input.GetAxisRaw("Mouse X") * mouseSensitivity;
            float my = Input.GetAxisRaw("Mouse Y") * mouseSensitivity;

            _yaw += mx;
            _pitch -= my;
            _pitch = Mathf.Clamp(_pitch, minPitch, maxPitch);

            transform.rotation = Quaternion.Euler(_pitch, _yaw, 0f);

            // W/S 调整相机高度（相对玩家站位的 Y 偏移），夹在 [minCameraHeight, maxCameraHeight]
            float heightDelta = 0f;
            if (Input.GetKey(KeyCode.W)) heightDelta += 1f;
            if (Input.GetKey(KeyCode.S)) heightDelta -= 1f;
            if (heightDelta != 0f)
            {
                cameraHeight = Mathf.Clamp(
                    cameraHeight + heightDelta * heightAdjustSpeed * Time.deltaTime,
                    minCameraHeight, maxCameraHeight);
            }

            // 跟随站位（允许玩家站位在 runtime 移动）
            var gm = BasketballGameManager.Instance;
            if (gm != null && gm.playerStand != null)
            {
                transform.position = gm.playerStand.position + Vector3.up * cameraHeight;
            }
        }

        public bool IsCursorLocked => _cursorLocked;

        private void LockCursor(bool locked)
        {
            _cursorLocked = locked;
            Cursor.lockState = locked ? CursorLockMode.Locked : CursorLockMode.None;
            Cursor.visible = !locked;
        }
    }
}
