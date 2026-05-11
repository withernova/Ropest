using UnityEngine;

namespace Ropest.Demo
{
    /// <summary>
    /// 极简的"伪 Cinemachine"相机控制器。
    /// 之所以不用 CinemachineVirtualCamera：
    ///   - 我们只需要"按段切换预设 + 段内时间归一化插值"，自己实现 5 行就够了；
    ///   - 不引入 prefab/asset 依赖，纯运行时构造，DemoDirector 在 Awake 里直接 new 即可。
    ///
    /// 用法：
    ///   1. DemoDirector 调 ApplyPreset(preset, duration) 切段时设定本段相机预设和持续时间；
    ///   2. 每帧 LateUpdate 自动按 ElapsedTime / Duration 做 SmoothStep 插值，更新主相机位姿；
    ///   3. 段间瞬切（不做 cross-fade，盲审视频里干净的硬切反而更"演示味"）。
    /// </summary>
    [DisallowMultipleComponent]
    public class DemoCameraRig : MonoBehaviour
    {
        public Camera TargetCamera;

        DemoCameraPreset _preset;
        float _duration;
        float _elapsed;
        bool _hasPreset;

        public void ApplyPreset(DemoCameraPreset preset, float duration)
        {
            _preset = preset;
            _duration = Mathf.Max(0.01f, duration);
            _elapsed = 0f;
            _hasPreset = true;

            if (TargetCamera != null)
                TargetCamera.fieldOfView = preset.FieldOfView;

            // 初始位姿立刻摆好，避免段切瞬间出现一帧"上一段位姿 + 新段 LookAt"的鬼影。
            UpdateCameraTransform(0f);
        }

        void LateUpdate()
        {
            if (!_hasPreset) return;
            _elapsed += Time.deltaTime;
            float t = Mathf.Clamp01(_elapsed / _duration);
            // SmoothStep：缓入缓出，比线性插值看起来更"电影感"。
            float s = t * t * (3f - 2f * t);
            UpdateCameraTransform(s);
        }

        void UpdateCameraTransform(float s)
        {
            if (TargetCamera == null) return;

            float distance = Mathf.Lerp(_preset.DistanceStart, _preset.DistanceEnd, s);
            float yawDeg = Mathf.Lerp(_preset.YawStart, _preset.YawEnd, s);
            float pitchDeg = _preset.Pitch;

            // 球坐标 → 笛卡尔
            float yaw = yawDeg * Mathf.Deg2Rad;
            float pitch = pitchDeg * Mathf.Deg2Rad;

            Vector3 dir = new Vector3(
                Mathf.Cos(pitch) * Mathf.Sin(yaw),
                Mathf.Sin(pitch),
                Mathf.Cos(pitch) * Mathf.Cos(yaw));

            Vector3 pos = _preset.LookAt + dir * distance;
            TargetCamera.transform.position = pos;
            TargetCamera.transform.rotation = Quaternion.LookRotation(_preset.LookAt - pos, Vector3.up);
        }
    }
}
