public class Singleton<T> where T : new()
{
    static T _instance;

    static Singleton()
    {
        _instance = new T();
    }

    public static T instance { 
        get { return _instance; } 
    }
}
