using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using Unity.Mathematics;
using UnityEngine;

public class Rope : MonoBehaviour
{
    [Header("环境参数")]
    [SerializeField] private Vector3 gravity = new(0, -9.8f, 0);
    public bool simulate = true;
    private bool ready = false;

    private Mesh ropeMesh;
    MeshCollider meshCollider;
    [Header("绳子参数")]
    [Range(0.1f, 100f)] public float length = 5;
    [Range(1, 100)] public int segments = 10;
    [Range(3, 30)] public int subdivision = 20;
    [Range(0.01f, 10f)] public float thickness = 1;
    public Vector3 meshOrigin;
    public Material ropeMaterial;
    
    // Solver
    public RopeXPBDSolver solver;

    [Header("游戏对象")]
    public PlayerController controlPoint;



    private void Awake()
    {
        ropeMesh = CreateRope(length, segments, subdivision, thickness, meshOrigin);
        if (simulate)
            StartCoroutine(InitSolver());
        gameObject.AddComponent<MeshFilter>();
        gameObject.AddComponent<MeshRenderer>();
    }

    IEnumerator InitSolver()
    {
        solver = new RopeXPBDSolver();
        solver.InitByMesh(ropeMesh, gravity, new RopeSolverInitData(subdivision, segments, transform.GetChild(1).gameObject.AddComponent<CapsuleCollider>()));
        controlPoint.ctrl = solver;
        ready = true;
        yield break;
    }

    private void Start()
    {
        GetComponent<MeshFilter>().mesh = ropeMesh;
        GetComponent<MeshRenderer>().material = ropeMaterial;
        meshCollider = gameObject.AddComponent<MeshCollider>();
        meshCollider.convex = true;
        meshCollider.isTrigger = true;
        meshCollider.sharedMesh = null;
        meshCollider.sharedMesh = ropeMesh;
        //gameObject.AddComponent<Rigidbody>().useGravity = false;
    }

    // 处理动作和更新mesh
    public void Update()
    {
        if (Input.GetMouseButtonDown(0))
        {
            Ray ray = Camera.main.ScreenPointToRay(Input.mousePosition);
            Physics.Raycast(ray, out RaycastHit hitInfo);
            Vector3 vertexPos = ray.origin + ray.direction * 2;
            solver.StartGrab(vertexPos, transform);
        }
        if(Input.GetMouseButtonUp(0))
        {
            Ray ray = Camera.main.ScreenPointToRay(Input.mousePosition);
            Physics.Raycast(ray, out RaycastHit hitInfo);
            if (solver.grabPoint != -1)
            {
                solver.EndGrab(hitInfo.point, transform);
            }
        }
        if (Input.GetMouseButton(0))
        {
            if (solver.grabPoint != -1)
            {
                Ray ray = Camera.main.ScreenPointToRay(Input.mousePosition);
                Physics.Raycast(ray, out RaycastHit hitInfo);
                Vector3 vertexPos = ray.origin + ray.direction * 2;
                solver.OnGrabbing(vertexPos, transform);
            }
        }


        if (ready)
        {            
            UpdateMesh();
        }
    }

    // 计算和模拟
    public void FixedUpdate()
    {
        if (!simulate || !ready)
            return;
        controlPoint.transform.localPosition = solver.pointPos[solver.ctrlIndex];
        // for (var index = 0; index < solver.pointPos.Length; index++)
        // {
        //     solver.direcN[index] = Vector3.zero;
        //     var localVertex = solver.pointPos[index];
        //     Vector3 worldVertex = transform.TransformPoint(localVertex);
        //     foreach (var collider in Physics.OverlapSphere(worldVertex, 0.05f))
        //     {
        //         if (collider.gameObject == gameObject || collider.isTrigger) continue;
        //
        //         Physics.Raycast(worldVertex, collider.ClosestPoint(worldVertex) - worldVertex, out RaycastHit info);
        //         solver.direcN[index] += info.normal;
        //     }
        // }
        solver.Simulate(Time.fixedDeltaTime);
        
    }

    private Mesh CreateRope(float len, int seg, int sub, float radius, Vector3 initPos)
    {
        Mesh mesh = new Mesh();
        mesh.name = "Rope";
        List<Vector3> vertices = new List<Vector3>();
        List<int> triangles = new List<int>();
        List<Vector2> uvs = new List<Vector2>();

        // 顶点
        for (int i = 0; i < seg + 1; i++)
        {
            Vector3 point = initPos + (len / seg) * new Vector3(i, 0, 0);
            float angleIncrement = 2 * Mathf.PI / sub;
            if (i == 0)
                vertices.Add(point);
            for (int j = 0; j < sub; j++)
            {
                float angle = j * angleIncrement;
                float x = point.x;
                float y = point.y + Mathf.Cos(angle) * radius;
                float z = point.z + Mathf.Sin(angle) * radius;
                vertices.Add(new Vector3(x, y, z));
            }
            if (i == seg)
                vertices.Add(point);
        }

        // 索引
        for (int i = 0; i < seg; i++)
        {
            if (i == 0)
            {
                for (int j = 0; j < sub; j++)
                {
                    int nextJ = (j + 1) % sub;
                    triangles.Add(0);
                    triangles.Add(1 + nextJ);
                    triangles.Add(1 + j);
                }
            }
            for (int j = 0; j < sub; j++)
            {
                int nextJ = (j + 1) % sub;
                triangles.Add(i * sub + 1 + j);
                triangles.Add(i * sub + 1 + nextJ);
                triangles.Add((i + 1) * sub + 1 + j);

                triangles.Add((i + 1) * sub + 1 + j);
                triangles.Add(i * sub + 1 + nextJ);
                triangles.Add((i + 1) * sub + 1 + nextJ);
            }
        }
        for (int j = 0; j < sub; j++)
        {
            int nextJ = (j + 1) % sub;
            triangles.Add(seg * sub + 1 + j);
            triangles.Add(seg * sub + 1 + nextJ);
            triangles.Add((seg + 1) * sub + 1);
        }


        // uv
        uvs.Add(new Vector2(0, 0));
        for (int i = 0; i <= seg; i++)
        {
            for (int j = 0; j < sub; j++)
            {
                float u, v;
                if (j <= sub / 2)
                {
                    u = j / sub * 2;
                }
                else
                {
                    u = (sub - j) / sub * 2;
                }
                v = i;
                uvs.Add(new Vector2(u, v));
            }
        }
        uvs.Add(new Vector2(0, 0));

        mesh.SetVertices(vertices);
        mesh.SetTriangles(triangles, 0);
        mesh.uv = uvs.ToArray();
        mesh.RecalculateNormals();
        mesh.RecalculateTangents();
        return mesh;
    }

    private void UpdateMesh()
    {
        solver.TryRendering(ropeMesh);
        //ropeMesh.vertices = solver.pos;
        ropeMesh.RecalculateBounds();
        ropeMesh.RecalculateNormals();

        meshCollider.sharedMesh = null;
        meshCollider.sharedMesh = ropeMesh;

    }

    private void OnDrawGizmos()
    {
        if(solver != null && solver.grabPoint != -1)
        {
            //Gizmos.DrawSphere(transform.TransformPoint(solver.pos[solver.grabPoint]), 1);
        }
    }
}
