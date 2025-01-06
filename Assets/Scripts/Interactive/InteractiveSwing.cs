using OutlineFx;
using System.Collections.Generic;
using UnityEngine;

public class InteractiveSwing : InteractiveBase
{
    private void Awake()
    {
        outline = gameObject.GetComponent<Outline>();
        outline.Color = Color.red;
        outline.enabled = false;
    }
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
