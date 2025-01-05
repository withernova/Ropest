using Cinemachine;
using System;
using System.Collections.Generic;
using TMPro;
using UnityEngine;
using UnityEngine.SceneManagement;
using UnityEngine.UI;

public class GameManager : MonoBehaviour
{
    public static GameManager Instance;

    public bool isUseTimer = true;
    public GameObject respawnPoint;
    public GameObject defeatedUI;
    private GameObject rope;
    public Transform deadPlane;

    public CinemachineFreeLook cam;
    public int score;
    public TextMeshProUGUI scoreTextUI;

    private void Awake()
    {
        Instance = this;
    }

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

        if(null != deadPlane)
        {
            deadPlane.transform.position = new Vector3(rope.transform.position.x, rope.transform.position.y, -5);
        }
    }

    public void OnGameOver()
    {
        defeatedUI.transform.DoPageAnimation(AnimationType.PumpOnce, true);
        scoreTextUI.text = score.ToString();
    }

    public void PlayAgain()
    {
        MapManager.Instance.previousMap.ForEach(a => ResourcesPool.Instance.ReturnOne(MapManager.Instance.maps[a.id], a.gameObject));
        MapManager.Instance.previousMap = new List<LevelMap>();
        SceneManager.LoadScene("Menu");
    }
}