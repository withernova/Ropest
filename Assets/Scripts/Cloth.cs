using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Serialization;
using Random = UnityEngine.Random;

public class Cloth : MonoBehaviour
{
    [Header("环境参数")]
    [SerializeField] private Vector3 gravity = new(0, -9.8f, 0);
    public bool simulate = true;
    private bool ready = false;

    private Mesh clothMesh;
    private MeshFilter clothMeshFilter;
    MeshCollider meshCollider;
    [Header("布料参数")]
    [Range(0.1f, 100f)] public float length = 5;
    [Range(0.1f, 100f)] public float width = 5;
    [Range(1, 100)] public int segments = 10;
    [Range(1, 100)] public int subdivision = 11;
    public Vector3 meshOrigin;
    public Material clothMaterial;
    
    //private List<Tuple<Vector3, Vector3>> vertice_force;
    //public List<Vector3> verticeInCollider;

    // Solver
    public ClothXPBDSolver solver;

    private void Awake()
    {
        clothMesh = CreateCloth(length, width, segments, subdivision, meshOrigin);
        if (simulate)
            StartCoroutine(InitSolver());
        gameObject.AddComponent<MeshFilter>();
        gameObject.AddComponent<MeshRenderer>();
    }

    IEnumerator InitSolver()
    {
        solver = new ClothXPBDSolver();
        solver.InitByMesh(clothMesh, gravity, new ClothSolverInitData(subdivision));
        ready = true;
        yield break;
    }

    private void Start()
    {
        clothMeshFilter = GetComponent<MeshFilter>();
        clothMeshFilter.mesh = clothMesh;
        GetComponent<MeshRenderer>().material = clothMaterial;
        meshCollider = gameObject.AddComponent<MeshCollider>();
        meshCollider.convex = true;
        meshCollider.isTrigger = true;
        meshCollider.sharedMesh = clothMesh;
        //gameObject.AddComponent<Rigidbody>().useGravity = false;
        //vertice_force = new List<Tuple<Vector3, Vector3>>();
    }

    // 处理动作和更新mesh
    public void Update()
    {
        if (Input.GetMouseButtonDown(0))
        {
            Ray ray = Camera.main.ScreenPointToRay(Input.mousePosition);
            Physics.Raycast(ray, out RaycastHit hitInfo);
            if(hitInfo.transform == transform)
            {
                solver.StartGrab(hitInfo.point, transform);
            }
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
            UpdateMesh();
        
        //vertice_force.Clear();
        //verticeInCollider.Clear();
    }

    // 计算和模拟
    public void FixedUpdate()
    {
        if (!simulate || !ready)
            return;
        for (var index = 0; index < clothMesh.vertices.Length; index++)
        {
            solver.direcN[index] = Vector3.zero;
            var localVertex = clothMesh.vertices[index];
            Vector3 worldVertex = clothMeshFilter.transform.TransformPoint(localVertex);
            foreach (var collider in Physics.OverlapSphere(worldVertex, 0.1f))
            {
                if (collider.gameObject == gameObject) continue;
                
                Physics.Raycast(worldVertex, collider.ClosestPoint(worldVertex) - worldVertex, out RaycastHit info);
                
                solver.direcN[index] += info.normal;
                if (collider.TryGetComponent<Rigidbody>(out var rb))
                {
                    Vector3 v1 = solver.vel[index];
                    Vector3 v2 = rb.linearVelocity;
                    float m1 = 1f / solver.invMass[index];
                    float m2 = rb.mass;
                    Vector3 v1_ = (v1 * (m1 - m2) + 2 * m2 * v2) / (m1 + m2);
                    Vector3 v2_ = (v2 * (m2 - m1) + 2 * m1 * v1) / (m1 + m2);
                    
                    //Debug.Log(m1);
                    //Debug.Log(v1);
                    //Debug.Log(v2_);
                    
                    
                    solver.vel[index] = v1_;
                    //solver.pos[index] = solver.prevPos[index] + solver.vel[index] * Time.fixedDeltaTime;
                    //pos[i] = prevPos[i] + vel[i] * dt;
                    rb.linearVelocity = v2_;
                }
            }
        }
        solver.Simulate(Time.fixedDeltaTime);
        
    }

    private Mesh CreateCloth(float len, float wid, int seg, int sub, Vector3 initPos)
    {
        Mesh mesh = new Mesh();
        mesh.name = "Cloth";
        List<Vector3> vertices = new List<Vector3>();
        List<int> triangles = new List<int>();
        List<Vector2> uvs = new List<Vector2>();

        // 顶点
        for (int i = 0; i < seg + 1; i++)
        {
            Vector3 point = initPos + (len / seg) * new Vector3(i, 0, 0);
            for (int j = 0; j < sub + 1; j++)
            {
                Vector3 point1 = point + (wid / sub) * new Vector3(0, 0, j);
                vertices.Add(point1);
            }
        }


        // 索引
        for (int i = 0; i < seg; i++)
        {
            for (int j = 0; j < sub; j++)
            {
                int nextJ = (j + 1) % (sub + 1);
                triangles.Add(i * (sub + 1) + j);
                triangles.Add(i * (sub + 1) + nextJ);
                triangles.Add((i + 1) * (sub + 1) + j);

                triangles.Add((i + 1) * (sub + 1) + j);
                triangles.Add(i * (sub + 1) + nextJ); 
                triangles.Add((i + 1) * (sub + 1) + nextJ);
            }
        }


        // uv
        for (int i = 0; i <= seg; i++)
        {
            for (int j = 0; j <= sub; j++)
            {
                float u, v;
                u = j;
                v = i;
                uvs.Add(new Vector2(u, v));
            }
        }

        mesh.SetVertices(vertices);
        mesh.SetTriangles(triangles, 0);
        mesh.uv = uvs.ToArray();
        mesh.RecalculateNormals();
        mesh.RecalculateTangents();
        return mesh;
    }
    
    private void UpdateMesh()
    {
        clothMesh.vertices = solver.pos;
        clothMesh.RecalculateNormals();
        clothMesh.RecalculateTangents();
        clothMesh.RecalculateBounds();

        meshCollider.sharedMesh = null;
        meshCollider.sharedMesh = clothMesh;

    }

    private void OnDrawGizmos()
    {
        if(solver != null && solver.grabPoint != -1)
        {
            //Gizmos.DrawSphere(transform.TransformPoint(solver.pos[solver.grabPoint]), 1);
        }
    }
}
