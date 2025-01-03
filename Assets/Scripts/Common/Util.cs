using UnityEngine;

public static class Util
{
    public static bool FloatEqual(float a, float b)
    {
        return Mathf.Abs(a - b) < 0.0001f;
    }
    public static bool DoubleEqual(double a, double b)
    {
        return Mathf.Abs((float)(a - b)) < 0.0001f;
    }
}
