using UnityEngine;

public interface IControllable
{
    public void SetMove(Vector3 move);

    public void Hang(Vector3 target);

    public void Grab(Vector3 target);
}
