using System;
using System.Collections.Generic;
using UnityEngine;

public class GameManager : SingletonForMonoBehaviour<GameManager>
{
    public bool isUseTimer = true;
    public GameObject respawnPoint;
    private GameObject rope;
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
    }

    private void Update()
    {
        if(isUseTimer)
            TimerManager.instance.Loop(Time.deltaTime);
    }

}
