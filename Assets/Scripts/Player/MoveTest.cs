using UnityEngine;

public class MoveTest : MonoBehaviour
{
    public float speed;
    // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Start()
    {
        GetComponent<Rigidbody>().linearVelocity = speed * Vector3.up;
    }

    // Update is called once per frame
    void Update()
    {
        
    }
}
