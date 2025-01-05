using Cinemachine;
using System;
using System.Collections.Generic;
using UnityEngine;

public class GameManager : SingletonForMonoBehaviour<GameManager>
{
    public bool isUseTimer = true;
    public GameObject respawnPoint;
    private GameObject rope;

    public CinemachineFreeLook cam;

    private void Start()
    {
        LoadRope();
    }

    public void LoadRope()
    {
        if(rope)
            Destroy(rope);
        rope = Instantiate(Resources.Load("Prefabs/Rope") as GameObject);
        rope.GetComponent<Rope>().LoadRope(respawnPoint.transform.position);
        cam.LookAt = rope.transform.GetChild(0);
        cam.Follow = rope.transform.GetChild(0);
    }

    private void Update()
    {
        if(isUseTimer)
            TimerManager.instance.Loop(Time.deltaTime);
    }

}