using UnityEngine;

namespace Ropest.Demo
{
    /// <summary>
    /// 单段相机预设。导演按段切换 DemoCameraRig，做 Cinemachine 风格的"环绕 + 拉远"运动。
    ///
    /// 运动定义（极坐标）：
    ///   - LookAt：视点（看向哪个世界点）。
    ///   - DistanceStart / DistanceEnd：随时间从起始距离平滑过渡到结束距离（拉远/拉近）。
    ///   - YawStart / YawEnd：水平环绕角度，单位度。从 YawStart 平滑过渡到 YawEnd。
    ///   - Pitch：俯仰角，单位度（正值为俯视）。
    ///   - FieldOfView：固定 FOV。
    /// </summary>
    [System.Serializable]
    public struct DemoCameraPreset
    {
        public Vector3 LookAt;
        public float DistanceStart;
        public float DistanceEnd;
        public float YawStart;
        public float YawEnd;
        public float Pitch;
        public float FieldOfView;

        /// <summary>静态镜头：起止距离/角度相同。</summary>
        public static DemoCameraPreset Static(Vector3 lookAt, float distance, float yaw, float pitch, float fov = 45f)
        {
            return new DemoCameraPreset
            {
                LookAt = lookAt,
                DistanceStart = distance,
                DistanceEnd = distance,
                YawStart = yaw,
                YawEnd = yaw,
                Pitch = pitch,
                FieldOfView = fov
            };
        }

        /// <summary>缓慢环绕：距离不变，yaw 从 a 到 b。</summary>
        public static DemoCameraPreset Orbit(Vector3 lookAt, float distance, float yawA, float yawB, float pitch, float fov = 45f)
        {
            return new DemoCameraPreset
            {
                LookAt = lookAt,
                DistanceStart = distance,
                DistanceEnd = distance,
                YawStart = yawA,
                YawEnd = yawB,
                Pitch = pitch,
                FieldOfView = fov
            };
        }

        /// <summary>边环绕边拉远。</summary>
        public static DemoCameraPreset OrbitAndZoom(Vector3 lookAt,
            float distA, float distB, float yawA, float yawB, float pitch, float fov = 45f)
        {
            return new DemoCameraPreset
            {
                LookAt = lookAt,
                DistanceStart = distA,
                DistanceEnd = distB,
                YawStart = yawA,
                YawEnd = yawB,
                Pitch = pitch,
                FieldOfView = fov
            };
        }

        /// <summary>反向预设：起点和终点互换。用于"可视化重播"阶段，让镜头从原终点反向走回起点。</summary>
        public DemoCameraPreset Reversed()
        {
            return new DemoCameraPreset
            {
                LookAt = LookAt,
                DistanceStart = DistanceEnd,
                DistanceEnd = DistanceStart,
                YawStart = YawEnd,
                YawEnd = YawStart,
                Pitch = Pitch,
                FieldOfView = FieldOfView
            };
        }
    }
}
