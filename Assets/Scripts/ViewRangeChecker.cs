using Cinemachine;
using UnityEngine;

public class ViewRangeChecker
{
    private Camera mainCamera;
    private CinemachineBrain cinemachineBrain;

    public ViewRangeChecker()
    {
        mainCamera = Camera.main;
        cinemachineBrain = mainCamera.GetComponent<CinemachineBrain>();
    }

    public bool IsInCameraView(Vector3 position, float threshold = 0.1f)
    {
        // 获取当前活动的虚拟相机
        ICinemachineCamera activeVCam = cinemachineBrain.ActiveVirtualCamera;

        // 转换为屏幕坐标
        Vector3 viewportPoint = mainCamera.WorldToViewportPoint(position);

        // 检查是否在视口范围内
        bool inView = viewportPoint.x >= -threshold && viewportPoint.x <= (1 + threshold) &&
                     viewportPoint.y >= -threshold && viewportPoint.y <= (1 + threshold) &&
                     viewportPoint.z > 0;

        return inView;
    }

    // 如果需要考虑距离和角度
    public bool IsInCameraRange(Vector3 position, float maxDistance = 10f, float maxAngle = 45f)
    {
        ICinemachineCamera activeVCam = cinemachineBrain.ActiveVirtualCamera;
        Transform vcamTransform = activeVCam.VirtualCameraGameObject.transform;

        // 检查距离
        float distance = Vector3.Distance(vcamTransform.position, position);
        if (distance > maxDistance) return false;

        // 检查角度
        Vector3 directionToTarget = (position - vcamTransform.position).normalized;
        float angle = Vector3.Angle(vcamTransform.forward, directionToTarget);
        if (angle > maxAngle) return false;

        return true;
    }

    // 完整的检查方法（包括遮挡检测）
    public bool IsFullyVisible(Vector3 position, float radius = 1f)
    {
        if (!IsInCameraView(position)) return false;

        // 检查遮挡
        Vector3 directionToCamera = mainCamera.transform.position - position;
        float distance = directionToCamera.magnitude;
        Ray ray = new Ray(position, directionToCamera);

        if (Physics.Raycast(ray, out RaycastHit hit, distance))
        {
            // 如果射线击中的不是相机，说明有遮挡
            if (hit.collider.gameObject != mainCamera.gameObject)
            {
                return false;
            }
        }

        return true;
    }

    // Debug绘制视野范围（在Scene视图中）
    void OnDrawGizmos()
    {
        if (cinemachineBrain == null || !Application.isPlaying) return;

        ICinemachineCamera activeVCam = cinemachineBrain.ActiveVirtualCamera;
        if (activeVCam == null) return;

        Transform vcamTransform = activeVCam.VirtualCameraGameObject.transform;

        // 绘制视锥体
        Gizmos.color = Color.yellow;
        float distance = 10f;
        float angle = 45f;

        Vector3 forward = vcamTransform.forward * distance;
        Vector3 right = Quaternion.Euler(0, angle, 0) * forward;
        Vector3 left = Quaternion.Euler(0, -angle, 0) * forward;

        Gizmos.DrawLine(vcamTransform.position, vcamTransform.position + forward);
        Gizmos.DrawLine(vcamTransform.position, vcamTransform.position + right);
        Gizmos.DrawLine(vcamTransform.position, vcamTransform.position + left);
    }
}