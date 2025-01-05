using System.Collections.Generic;
using UnityEngine;

public class InteractiveSwing : InteractiveBase
{
    public override void OnStart()
    {
        activating = true;
    }

    public override void OnUpdate(float deltaTime)
    {
    }

    public override void EndInteractive()
    {
        activating = false;
    }
}
