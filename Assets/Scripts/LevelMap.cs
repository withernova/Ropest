using DG.Tweening;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UnityEngine;


public class LevelMap : MonoBehaviour
{
    public int time;
    public int id;
    public Timer dropMapCountDown;
    public Timer onDropMap;
    public bool droped;
    public int hardLevel;

    public virtual void Init(int i)
    {
        id = i;
        time = hardLevel * 20;
        dropMapCountDown?.Stop();
        dropMapCountDown = TimerManager.instance.CreateTimer(time, 1, () => MapManager.Instance.DropMap(this));
        droped = false;
        transform.position = new Vector3(transform.position.x, -1);

        Material mat = MapManager.Instance.materials[Random.Range(0, MapManager.Instance.materials.Count)];
        transform.Find("object").GetComponentsInChildren<MeshRenderer>().ToList().ForEach(x => x.material = mat);

        //GetComponentsInChildren<InteractiveBase>().ToList().ForEach(inter => inter.EndInteractive());
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
        if (other.gameObject.CompareTag("Player") && other.transform.position.x < GetLasPos().x)
        {
            //dropMapCountDown.LastTime = Mathf.Min(dropMapCountDown.LastTime, 1f);
            dropMapCountDown.Stop();
            MapManager.Instance.DropMap(this);
            other.transform.parent.GetComponent<Rope>().solver.ReleaseAll();
            GameManager.Instance.score += 1;
        }
    }
}