using UnityEngine;

public class MoveTest : MonoBehaviour
{
    public float speed;
    void Start()
    {
        GetComponent<Rigidbody>().linearVelocity = speed * Vector3.up;
    }

    void Update()
    {
        
    }
}
