using UnityEngine;

public interface IControllable
{
    public void SetMove(Vector3 move);

    public void Hang(Vector3 controllable);

    public void Grab();
}
