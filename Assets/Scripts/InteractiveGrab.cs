using UnityEngine;

public class InteractiveGrab : InteractiveBase
{
    RigidbodyXPBD rigidbodyXPBD;

    private void Start()
    {
        rigidbodyXPBD = GetComponent<RigidbodyXPBD>();
    }

    public float GetMass()
    {
        return rigidbodyXPBD.mass;
    }

    public void SetV(Vector3 v)
    {
        rigidbodyXPBD.linearVelocity = v;
    }

    public override void OnStart()
    {
        GetComponent<Rigidbody>().isKinematic = true;
    }

    public override void OnUpdate(float deltaTime)
    {

    }

    public override void EndInteractive()
    {
        GetComponent<Rigidbody>().isKinematic = false;
        rigidbodyXPBD.rb.linearVelocity = rigidbodyXPBD.linearVelocity;
    }
}
