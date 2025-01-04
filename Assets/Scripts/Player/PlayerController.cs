using System;
using Cinemachine;
using NUnit.Framework;
using UnityEngine;

public class PlayerController : MonoBehaviour
{ 
    public CinemachineFreeLook freeLookCamera;
    private SphereCollider groundCheckCollider;
    private float speed = 5f;
    private float jumpSpeed = 1f;
    private Vector3 move;
    public IControllable ctrl;

    public InteractiveCaster caster;
    private void Start()
    {
        freeLookCamera.Follow = transform;
        freeLookCamera.LookAt = transform;
        caster = GetComponent<InteractiveCaster>();

        //groundCheckCollider = GetComponent<SphereCollider>();
    }

    private void Update()
    {
        //// ×¥È¡µÈ
        if (ctrl != null)
        {
            //Transform transform;

            if (Input.GetKeyDown(KeyCode.E))
            {
                InteractiveGrab interactive = caster.TriggerInteractiveUndeploy<InteractiveGrab>();
                if (interactive != null)
                {
                    Debug.Log("a");
                    ctrl.Grab(interactive);
                }
            }

            if (Input.GetKeyDown(KeyCode.Q))
            {
                InteractiveSwing interactive = caster.TriggerInteractiveUndeploy<InteractiveSwing>();
                Debug.Log(interactive);
                if (interactive != null)
                    ctrl.Swing(interactive);
            }
        }
    }

    private void FixedUpdate()
    {
        move = Vector3.zero;
        float horizontal = 0f;
        float vertical = 0f;
        horizontal = Input.GetAxis("Horizontal");
        vertical = Input.GetAxis("Vertical");
        Transform cameraTransform = freeLookCamera.VirtualCameraGameObject.transform;
        Vector3 forward = Camera.main.transform.forward;
        Vector3 horizontalForward = new Vector3(forward.x, 0, forward.z).normalized;
        Vector3 right = Vector3.Cross(horizontalForward, Vector3.up);
        //move += Time.deltaTime * new Vector3(1, 0, 0);
        //move += Time.deltaTime * speed * vertical * horizontalForward + Time.deltaTime * speed * horizontal * right;
        move += Time.deltaTime * speed * vertical * horizontalForward;
        // if (groundCheckCollider.CompareTag("Ground"))
        // {s
        //     if (Input.GetKeyDown(KeyCode.Space))
        //     {
        //         move += jumpSpeed * Vector3.up;
        //     }
        // }



        if (ctrl != null)
        {
            ctrl.SetMove(move);
        }
    }


}
