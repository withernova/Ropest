using System;
using Cinemachine;
using NUnit.Framework;
using Unity.VisualScripting;
using UnityEngine;

public class PlayerController : MonoBehaviour
{
    private SphereCollider groundCheckCollider;
    private float speed = 1f;
    private float jumpSpeed = 1f;
    private Vector3 move;
    public IControllable ctrl;

    public InteractiveCaster caster;

    private void Start()
    {
        caster = GetComponent<InteractiveCaster>();

        caster.condition += item => !item.activating;
        //groundCheckCollider = GetComponent<SphereCollider>();
    }
    InteractiveGrab showGrab;
    InteractiveSwing showSwing;


    private void Update()
    {
        InteractiveGrab show1 = caster.TriggerInteractiveUndeploy<InteractiveGrab>();
        if(show1 != showGrab) {
            if (showGrab != null) showGrab.outline.enabled = false;
            showGrab = show1;
            if (show1 != null) show1.outline.enabled = true;
        }
        InteractiveSwing show2 = caster.TriggerInteractiveUndeploy<InteractiveSwing>();
        if (showSwing != show2)
        {
            if(showSwing != null) showSwing.outline.enabled = false;
            showSwing = show2;
            if (show2 != null) show2.outline.enabled = true;
        }

        if (ctrl != null)
        {
            //Transform transform;

            if (Input.GetKeyDown(KeyCode.E))
            {
                //InteractiveGrab interactive = caster.TriggerInteractiveUndeploy<InteractiveGrab>();
                //Debug.Log(interactive);

                ctrl.Grab(showGrab);

            }

            if (Input.GetKeyDown(KeyCode.Q))
            {
                //InteractiveSwing interactive = caster.TriggerInteractiveUndeploy<InteractiveSwing>();
                ctrl.Swing(showSwing);
            }

            if (Input.GetKeyDown(KeyCode.R))
            {
                ctrl.ReleaseAll();
            }

            if (Input.GetKeyDown(KeyCode.Alpha1))
            {
                ctrl.SwitchCtrl(1);
            }
            if (Input.GetKeyDown(KeyCode.Alpha3))
            {
                ctrl.SwitchCtrl(-1);
            }
            if (Input.GetKey(KeyCode.Space))
            {
                ctrl.Stop();
            }
        }

        if(transform.position.y < -5)
        {
            //Debug.Log(transform.position.y);
            GameManager.Instance.OnGameOver();
        }
    }

    private void FixedUpdate()
    {
        move = Vector3.zero;
        float horizontal = 0f;
        float vertical = 0f;
        horizontal = Input.GetAxis("Horizontal");
        vertical = Input.GetAxis("Vertical");
        Vector3 forward = Camera.main.transform.forward;
        Vector3 horizontalForward = new Vector3(forward.x, 0, forward.z).normalized;
        Vector3 right = Vector3.Cross(Vector3.up, horizontalForward);
        //move += Time.deltaTime * new Vector3(1, 0, 0);
        //move += Time.deltaTime * speed * vertical * horizontalForward + Time.deltaTime * speed * horizontal * right;
        move += Time.deltaTime * speed * (vertical * horizontalForward + horizontal * right).normalized;
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
