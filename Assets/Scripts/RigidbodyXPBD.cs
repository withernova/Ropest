using System;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Serialization;

public class RigidbodyXPBD : MonoBehaviour
{
    public Vector3 linearVelocity;
    private Rigidbody rb;
    private void Start()
    {
        rb = GetComponent<Rigidbody>();
        linearVelocity = Vector3.zero;
    }
    void Update()
    {
        rb.linearVelocity = linearVelocity;
    }
}
