using DG.Tweening;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using UnityEngine;


public class LevelMap : MonoBehaviour
{
    public int time;
    public int id;
    public Timer dropMapCountDown;
    public bool droped;

    public void Init(int i, int lastTime)
    {
        id = i;
        time = lastTime;
        dropMapCountDown = TimerManager.instance.CreateTimer(lastTime, 1, () => MapManager.Instance.DropMap(this));
        droped = false;
        transform.position = new Vector3(transform.position.x, 0);
    }

    public Vector3 GetLasPos()
    {
        return this.transform.Find("endingPos").position;
    }

    public void SetBeginningPos(float pos)
    {
        Vector3 newPos = transform.position;
        newPos.x = pos + transform.position.x - this.transform.Find("beginningPos").position.x;
        transform.position = newPos;
        //float x = this.transform.Find("beginningPos").position.x;
        //this.transform.Find("beginningPos").position = new Vector3(pos, 0, 0);
        //this.transform.Find("object").position += new Vector3(pos - x, 0, 0);
    }

    private void OnTriggerEnter(Collider other)
    {
        if (other.gameObject.CompareTag("Player"))
        {
            if (!dropMapCountDown.IsRunning) dropMapCountDown.Start();
        }
    }

    private void OnTriggerExit(Collider other)
    {
        if (other.gameObject.CompareTag("Player") && !droped)
        {
            //dropMapCountDown.LastTime = Mathf.Min(dropMapCountDown.LastTime, 1f);
            dropMapCountDown.Stop();
            MapManager.Instance.DropMap(this);
        }
    }
}