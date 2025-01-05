using UnityEngine;
using UnityEngine.SceneManagement;

public class MenuManager : MonoBehaviour
{
    public Transform targetmodel;
    private GameObject[] tips;//需要旋转的物体
 
    private void Awake()
    {
        tips = GameObject.FindGameObjectsWithTag($"tip");
    }
 
    private void Update()
    {
        foreach(GameObject tip in tips)
        {
            tip.transform.eulerAngles = new Vector3(targetmodel.eulerAngles.x, targetmodel.eulerAngles.y, 0);
        }
    }

    public void StartGame()
    {
        SceneManager.LoadScene("KUMO");
    }

    public void Menu()
    {
        SceneManager.LoadScene("Menu");
    }
}
