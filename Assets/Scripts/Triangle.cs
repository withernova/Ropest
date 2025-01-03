using System.Collections.Generic;
using System.Linq;
using TMPro;
using Unity.Mathematics;
using Unity.VisualScripting;
using UnityEngine;

public class Triangle
{
    public int[] pointsId;
    public float volume;
    public float lambda = 0;
    public Matrix4x4 initialMatrix; // 存储初始形状矩阵 点的坐标 质心 初始
    public Vector3 centroid;

    public Triangle(int[] pointsId)
    {
        this.pointsId = pointsId;
    }
   
}