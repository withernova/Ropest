using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;

public class UnityMeshApplier : MonoBehaviour
{
    [Header("环境参数")]
    [SerializeField] private Vector3 gravity = new(0, -9.8f, 0);
    public bool simulate = true;
    private bool ready = false;

    private Mesh ropeMesh;
    public Material ropeMaterial;


    // Solver
    private ClothXPBDSolver solver;

    private void Awake()
    {
        ropeMesh = GetComponent<MeshFilter>().mesh;
        if (simulate)
            StartCoroutine(InitSolver());
    }

    IEnumerator InitSolver()
    {
        solver = new ClothXPBDSolver();
        //solver.Init(ropeMesh, gravity);
        ready = true;
        yield break;
    }

    private void Start()
    {
        GetComponent<MeshRenderer>().material = ropeMaterial;
    }

    // 处理动作和更新mesh
    public void Update()
    {

        if (simulate && ready)
            UpdateMesh();
    }

    // 计算和模拟
    public void FixedUpdate()
    {
        if (!simulate || !ready)
            return;
        solver.Simulate(Time.fixedDeltaTime);
    }

    private void UpdateMesh()
    {
        ropeMesh.vertices = solver.pos;
        ropeMesh.RecalculateBounds();
        ropeMesh.RecalculateNormals();
    }
}