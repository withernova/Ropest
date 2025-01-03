using System;

//管理所有定时器
public class TimerManager : Singleton<TimerManager>
{
    public  event Action<float> TimerLoopCallback;

    public Timer CreateTimer(float deltaTime, int repeatTimes, Action callback)
    {
        Timer timer = new Timer();
        timer.DeltaTime = deltaTime;
        timer.RepeatTime = repeatTimes;
        timer.Callback = callback;

        return timer;
    }
    public Timer CreateAndStartTimer(float deltaTime, int repeatTimes, Action callback)
    {
        Timer timer = new Timer();
        timer.DeltaTime = deltaTime;
        timer.RepeatTime = repeatTimes;
        timer.Callback = callback;
        timer.Start();

        return timer;
    }

    public void Loop(float deltaTime)
    {
        if(TimerLoopCallback != null)
        {
            TimerLoopCallback(deltaTime);
        }
    }

}

public class Timer
{
    public float DeltaTime;
    public float LastTime
    {
        get => DeltaTime - _duringTime;
        set => _duringTime = DeltaTime - value;
    }
    public int RepeatTime;
    public Action Callback;

    float _duringTime;
    int _repeatedTimes;

    public bool IsRunning = false;

    private void reset()
    {
        _duringTime = 0;
        _repeatedTimes = 0;
        
        Pause();
    }

    public void Start()
    {
        reset();
        IsRunning = true;
        TimerManager.instance.TimerLoopCallback += Loop;
    }

    public void Pause()
    {
        IsRunning = false;
        TimerManager.instance.TimerLoopCallback -= Loop;
    }

    public void Stop()
    {
        Pause();
        reset();
    }

    public void Loop(float deltaTime)
    {
        _duringTime += deltaTime;
        if(_duringTime > DeltaTime || Util.FloatEqual(_duringTime, DeltaTime))
        {
            ++_repeatedTimes;
            _duringTime -= DeltaTime;

            if(Callback != null)
            {
                Callback();
            }

            if(RepeatTime > 0 && _repeatedTimes >= RepeatTime)
            {
                Stop();
            }
        }
    }
}

