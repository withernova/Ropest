using NUnit.Framework;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;

public class PlayerControllerByForce : MonoBehaviour
{
    public Transform curTrans; // 当前操纵小节的位置 最好在rope写个接口好点儿
    public int curSectionIndex = 0;// 当前操纵小节的下标
    public float detectRaius = 10f;
    public Camera mainCamera;


    public Rope rope;

    private void Start()
    {
        mainCamera = Camera.main;
        rope = GetComponent<Rope>();
    }


    void Update()
    {
        Movement();
    }

    private void Movement()
    {
        Vector3 force = new Vector3();
        float horizontal = 0f;
        float vertical = 0f;
        horizontal = Input.GetAxis("Horizontal");
        vertical = Input.GetAxis("Vertical");
        //force += new Vector3()；
        Vector3 front = mainCamera.transform.forward;
        Vector3 right = mainCamera.transform.right;
        front.y = 0f;
        right.y = 0f;
        front.Normalize();
        right.Normalize();
        force = front * vertical + right * horizontal;
        //Debug.Log($"这就是力量{force}");
        rope.solver.AddForce(curSectionIndex,force * 10);

        var hangables = DetectInteractItem<Hangabletem>(rope.transform.TransformPoint(rope.solver.pointPos[0]));
        var checker = new ViewRangeChecker();
        Hangabletem curHangable = null;
        for(int i = 0;i <hangables.Count;i++)
        {
            //Debug.Log(checker.IsInCameraView(hangables[i].Key.transform.position));
            if (checker.IsInCameraView(hangables[i].Key.transform.position))
            {
                curHangable = hangables[i].Key;
                break;
            }
        }

        if (Input.GetKeyDown(KeyCode.E))
        {
            if (curHangable != null)
            {
                StartHanging(curHangable, rope.transform.TransformPoint(rope.solver.pointPos[0]));
            }
        }

        if (Input.GetKeyUp(KeyCode.E))
        {
            StopHanging();
            //EndGrab();
        }
    }



    public void StopHanging()
    {

    }

    public void StartHanging(Hangabletem item,Vector3 oriPos)
    {
        Vector3 target = item.GetTargetPos(oriPos);

        // TODO：总之做个移动动画吧 下面的协程里面写
        StartCoroutine(MoveTowards(target,1f));
        
    }

    public IEnumerator MoveTowards(Vector3 tar,float duration)
    {
        Debug.Log("catch");

        while (Input.GetKey(KeyCode.E))
        {
            rope.solver.AddForce(0, (rope.transform.InverseTransformPoint(tar) - rope.solver.pointPos[0]) / (Time.fixedDeltaTime * Time.fixedDeltaTime) * 0.8f);
            Debug.Log(rope.transform.InverseTransformPoint(tar) - rope.solver.pointPos[0]);
            yield return new WaitForFixedUpdate();
        }
    }

    public List<KeyValuePair<T, float>> DetectInteractItem<T>(Vector3 oriPos) where T : Component
    {
        //返回检测是否碰得到可抓获的物体 通过layermask

        List<KeyValuePair<T, float>> items = new List<KeyValuePair<T, float>>();

        foreach (var collider in Physics.OverlapSphere(oriPos, detectRaius))
        {
            // TODO 返回时需要在屏幕上给出提示 并根据屏幕方向 显示最近且可选中的目标 构建在列表内

            var component = collider.GetComponent<T>();
            if (component != null)
            {
                float distance = (oriPos - collider.ClosestPoint(oriPos)).magnitude;
                items.Add(new(component, distance));
            }
        }
        var sortedList = items.OrderByDescending(item => item.Value).ToList();

        return sortedList;
    }
}
