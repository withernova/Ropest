using DG.Tweening;
using UnityEngine;

public class InteractiveKite : InteractiveBase
{
    public Transform platform;
    bool isTriggered = false;
    public override void OnStart()
    {
        if (isTriggered) return;
        isTriggered = true;
        platform.DOLocalMoveY(-0.098f, 2);
        transform.parent.GetComponent<Rigidbody>().isKinematic = true;
    }

    public override void OnUpdate(float deltaTime)
    {

    }

    public override void EndInteractive()
    {
        platform.position = new Vector3(platform.position.x, 3f);
        transform.position = platform.position + new Vector3(2.5f, 0, 0);
        transform.parent.GetComponent<Rigidbody>().isKinematic = false;
    }
    // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Start()
    {

    }

    private void Update()
    {
        if (transform.position.y < -4)
        {
            transform.parent.position = platform.position + new Vector3(2f, 0, 0);
            transform.parent.GetComponent<Rigidbody>().linearVelocity = Vector3.zero;
        }
    }
}
