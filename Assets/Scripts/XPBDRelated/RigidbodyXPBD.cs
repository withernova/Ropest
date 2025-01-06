using System;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Serialization;

[RequireComponent(typeof(Rigidbody))]
public class RigidbodyXPBD : MonoBehaviour
{
    public Vector3 linearVelocity;
    public Rigidbody rb;
    public float mass = 1f;

    private void Awake()
    {
        rb = GetComponent<Rigidbody>();
        linearVelocity = Vector3.zero;

    }

    private void Start()
    {
    }
    void Update()
    {
        //rb.linearVelocity = linearVelocity;
    }
}
