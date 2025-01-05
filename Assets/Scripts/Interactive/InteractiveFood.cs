using DG.Tweening;
using UnityEngine;

public class InteractiveFood : InteractiveBase
{
    bool isTriggered;
    public Transform blockingWall;

    public override void OnStart()
    {
        if (isTriggered) return; 
        isTriggered = true;
        blockingWall.DOLocalMoveY(-0.89f, 2);
    }

    public override void OnUpdate(float deltaTime)
    {

    }

    public override void EndInteractive()
    {
        blockingWall.transform.position = new Vector3(blockingWall.transform.position.x, 1);
    }
}
