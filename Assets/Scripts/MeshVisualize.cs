using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class MeshVisualize : MonoBehaviour
{
    private XPBDSolver solver;
    private Rope ropeComponent;
    private Color[] triangleColors; // �洢ÿ�������ε���ɫ

    public float lineWidth = 0.01f;
    public bool showTriangles = true;
    public bool showEdges = true;
    public float alpha = 0.6f; // ͸���ȿ���

    void Start()
    {
        ropeComponent = GetComponent<Rope>();
        StartCoroutine(WaitForSolver());
    }

    private IEnumerator WaitForSolver()
    {
        while (ropeComponent.solver == null)
        {
            yield return null;
        }
        solver = ropeComponent.solver;
        InitializeTriangleColors();
    }

    private void InitializeTriangleColors()
    {
        if (solver == null || solver.triangles == null) return;

        triangleColors = new Color[solver.triangles.Length];
        for (int i = 0; i < solver.triangles.Length; i++)
        {
            // ʹ��HSV��ɫ�ռ����ɸ����ȵ���ɫ�ֲ�
            float hue = (float)i / solver.triangles.Length; // 0��1֮���ֵ
            triangleColors[i] = Color.HSVToRGB(hue, 0.8f, 0.8f);
            triangleColors[i].a = alpha; // ����͸����
        }
    }

    void Update()
    {
        if (solver == null || solver.triangles == null || solver.pos == null || triangleColors == null) return;

        HashSet<string> drawnEdges = new HashSet<string>();

        for (int i = 0; i < solver.triangles.Length; i++)
        {
            Triangle triangle = solver.triangles[i];
            Vector3 p1 = transform.TransformPoint(solver.pos[triangle.pointsId[0]]);
            Vector3 p2 = transform.TransformPoint(solver.pos[triangle.pointsId[1]]);
            Vector3 p3 = transform.TransformPoint(solver.pos[triangle.pointsId[2]]);

            // ʹ��Debug.DrawLine���������α�Ե
            if (showEdges)
            {
                Color edgeColor = triangleColors[i];
                edgeColor.a = 1f; // ��Եʹ�ò�͸��ɫ

                string edgeKey = GetEdgeKey(p1, p2);
                if (!drawnEdges.Contains(edgeKey))
                {
                    Debug.DrawLine(p1, p2, edgeColor);
                    drawnEdges.Add(edgeKey);
                }

                edgeKey = GetEdgeKey(p2, p3);
                if (!drawnEdges.Contains(edgeKey))
                {
                    Debug.DrawLine(p2, p3, edgeColor);
                    drawnEdges.Add(edgeKey);
                }

                edgeKey = GetEdgeKey(p3, p1);
                if (!drawnEdges.Contains(edgeKey))
                {
                    Debug.DrawLine(p3, p1, edgeColor);
                    drawnEdges.Add(edgeKey);
                }
            }
        }
    }

    void OnDrawGizmos()
    {
        if (!Application.isPlaying || solver == null || solver.triangles == null ||
            solver.pos == null || triangleColors == null) return;

        HashSet<string> drawnEdges = new HashSet<string>();

        for (int i = 0; i < solver.triangles.Length; i++)
        {
            Triangle triangle = solver.triangles[i];
            Vector3 p1 = transform.TransformPoint(solver.pos[triangle.pointsId[0]]);
            Vector3 p2 = transform.TransformPoint(solver.pos[triangle.pointsId[1]]);
            Vector3 p3 = transform.TransformPoint(solver.pos[triangle.pointsId[2]]);

            if (showTriangles)
            {
                Gizmos.color = triangleColors[i];
                DrawTriangleFace(p1, p2, p3);
            }

            if (showEdges)
            {
                Color edgeColor = triangleColors[i];
                edgeColor.a = 1f; // ��Եʹ�ò�͸��ɫ
                Gizmos.color = edgeColor;

                DrawEdgeIfNotDrawn(p1, p2, drawnEdges);
                DrawEdgeIfNotDrawn(p2, p3, drawnEdges);
                DrawEdgeIfNotDrawn(p3, p1, drawnEdges);
            }
        }
    }

    private void DrawEdgeIfNotDrawn(Vector3 p1, Vector3 p2, HashSet<string> drawnEdges)
    {
        string edgeKey = GetEdgeKey(p1, p2);
        if (!drawnEdges.Contains(edgeKey))
        {
            Gizmos.DrawLine(p1, p2);
            drawnEdges.Add(edgeKey);
        }
    }

    private string GetEdgeKey(Vector3 p1, Vector3 p2)
    {
        return p1.x < p2.x ?
            $"{p1.x},{p1.y},{p1.z}-{p2.x},{p2.y},{p2.z}" :
            $"{p2.x},{p2.y},{p2.z}-{p1.x},{p1.y},{p1.z}";
    }

    private void DrawTriangleFace(Vector3 p1, Vector3 p2, Vector3 p3)
    {
        // ������������
        Vector3 center = (p1 + p2 + p3) / 3f;
        Vector3[] vertices = new Vector3[] { p1, p2, p3 };

        // ���ƴ����ĵ�������������������Ч��
        for (int i = 0; i < 3; i++)
        {
            Gizmos.DrawLine(center, vertices[i]);
        }

        // ���Ʊ߿�
        Gizmos.DrawLine(p1, p2);
        Gizmos.DrawLine(p2, p3);
        Gizmos.DrawLine(p3, p1);
    }

    // �ṩ����������ɫ�ķ���
    public void RegenerateColors()
    {
        InitializeTriangleColors();
    }
}