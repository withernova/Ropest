using Cinemachine;
using UnityEngine;

public class PlayerController : MonoBehaviour
{ 
    public CinemachineFreeLook freeLookCamera;
    private SphereCollider groundCheckCollider;
    private float speed;
    private float jumpSpeed;
    private Vector3 move;
    public IControllable ctrl;
    private void Start()
    {
        freeLookCamera.Follow = transform;
        freeLookCamera.LookAt = transform;
        groundCheckCollider = GetComponent<SphereCollider>();
    }

    private void Update()
    {
        move = Vector3.zero;
        Transform cameraTransform = freeLookCamera.VirtualCameraGameObject.transform;
        Vector3 forward = cameraTransform.forward;
        Vector3 horizontalDirection = new Vector3(forward.x, 0, forward.z).normalized;
        move += speed * horizontalDirection;
        if (groundCheckCollider.CompareTag("Ground"))
        {
            if (Input.GetKeyDown(KeyCode.Space))
            {
                move += jumpSpeed * Vector3.up;
            }
        }

        if (ctrl != null)
        {
            ctrl.SetMove(move);
        }
    }
}
