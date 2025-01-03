using System;
using System.Collections.Generic;
using UnityEngine;

public class GameManager : SingletonForMonoBehaviour<GameManager>
{

    private void Update()
    {
        TimerManager.instance.Loop(Time.deltaTime);
    }

}
