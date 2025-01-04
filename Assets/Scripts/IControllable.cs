using UnityEngine;

public interface IControllable
{
    public void SetMove(Vector3 move);

    public void Swing(InteractiveSwing target);

    public void Grab(InteractiveGrab target);

    public void LoseControl();
}
