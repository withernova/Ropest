using UnityEngine;

public class InteractiveFood : InteractiveBase
{
    bool isTriggered;

    public override void OnStart()
    {
        if (isTriggered) return; 
        isTriggered = true;
        Debug.Log("обр╩╧ь");
    }

    public override void OnUpdate(float deltaTime)
    {

    }
}
