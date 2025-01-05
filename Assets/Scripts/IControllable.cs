using UnityEngine;

public interface IControllable
{
    public void SetMove(Vector3 move);

    public void Swing(InteractiveSwing target);

    public void Grab(InteractiveGrab target);

    public void SwitchCtrl(int ctrlIndex);

    public void ReleaseAll();
    public void LoseControl();
}
