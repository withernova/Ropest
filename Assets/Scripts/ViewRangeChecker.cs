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
        ICinemachineCamera activeVCam = cinemachineBrain.ActiveVirtualCamera;
        Vector3 viewportPoint = mainCamera.WorldToViewportPoint(position);

        bool inView = viewportPoint.x >= -threshold && viewportPoint.x <= (1 + threshold) &&
                     viewportPoint.y >= -threshold && viewportPoint.y <= (1 + threshold) &&
                     viewportPoint.z > 0;

        return inView;
    }

    public bool IsInCameraRange(Vector3 position, float maxDistance = 10f, float maxAngle = 45f)
    {
        ICinemachineCamera activeVCam = cinemachineBrain.ActiveVirtualCamera;
        Transform vcamTransform = activeVCam.VirtualCameraGameObject.transform;

        float distance = Vector3.Distance(vcamTransform.position, position);
        if (distance > maxDistance) return false;

        Vector3 directionToTarget = (position - vcamTransform.position).normalized;
        float angle = Vector3.Angle(vcamTransform.forward, directionToTarget);
        if (angle > maxAngle) return false;

        return true;
    }

    public bool IsFullyVisible(Vector3 position, float radius = 1f)
    {
        if (!IsInCameraView(position)) return false;

        Vector3 directionToCamera = mainCamera.transform.position - position;
        float distance = directionToCamera.magnitude;
        Ray ray = new Ray(position, directionToCamera);

        if (Physics.Raycast(ray, out RaycastHit hit, distance))
        {
            if (hit.collider.gameObject != mainCamera.gameObject)
            {
                return false;
            }
        }

        return true;
    }

    void OnDrawGizmos()
    {
        if (cinemachineBrain == null || !Application.isPlaying) return;

        ICinemachineCamera activeVCam = cinemachineBrain.ActiveVirtualCamera;
        if (activeVCam == null) return;

        Transform vcamTransform = activeVCam.VirtualCameraGameObject.transform;

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
