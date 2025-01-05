using System;
using System.Collections.Generic;
using System.Linq;
using Cinemachine;
using NUnit.Framework;
using UnityEngine;

public class PointData
{
    public bool isActivating = false;
    public InteractiveBase interactiveItem;
    public float originMass;
    public int index;

    public PointData(float oriM,int i)
    {
        originMass = oriM;
        index = i;  
        PointData.datas.Add(this);
    }
   
    public bool SetActive(InteractiveBase baseI)
    {
        if(isActivating)
        {
            return false;
        }
        interactiveItem = baseI;
        baseI.OnStart();
        isActivating = true;
        return true;
    }
    public float Reset()
    {
        this.interactiveItem.EndInteractive();
        this.interactiveItem = null;
        this.isActivating = false;
        return originMass;
    }

    public static List<PointData> datas = new List<PointData>();

    public static List<PointData> GetAllActivated()
    {
        return datas.Where(data => data.isActivating).ToList();
    }
}