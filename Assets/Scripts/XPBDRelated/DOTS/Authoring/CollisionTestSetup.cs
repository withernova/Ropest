using UnityEngine;

public class CollisionTestSetup : MonoBehaviour
{
    [Header("球体参数")]
    public float sphereRadius = 0.5f;
    public Vector3 sphereStartPos = new Vector3(0, -2f, 0);
    public Vector3 sphereEndPos = new Vector3(3f, -2f, 0);
    public float sphereMoveSpeed = 1.5f;

    [Header("方体参数")]
    public Vector3 boxSize = new Vector3(1f, 1f, 1f);
    public Vector3 boxStartPos = new Vector3(-2f, -3f, 0);
    public Vector3 boxEndPos = new Vector3(2f, -3f, 0);
    public float boxMoveSpeed = 1.0f;

    [Header("运动模式")]
    [Tooltip("运动时间（秒）")]
    public float moveDuration = 2.0f;
    [Tooltip("停顿时间（秒）")]
    public float pauseDuration = 1.0f;

    [Header("材质")]
    public Material sphereMaterial;
    public Material boxMaterial;

    // 运行时对象
    private GameObject _sphereObj;
    private GameObject _boxObj;

    // 运动状态
    private float _sphereTimer;
    private float _boxTimer;
    private bool _sphereMovingForward = true;
    private bool _boxMovingForward = true;
    private bool _spherePaused = false;
    private bool _boxPaused = false;
    private float _spherePauseTimer;
    private float _boxPauseTimer;

    void Start()
    {
        CreateSphere();
        CreateBox();

        // 确保场景中有碰撞管理器
        if (AnalyticalColliderManager.Instance == null)
        {
            var managerObj = new GameObject("AnalyticalColliderManager");
            managerObj.AddComponent<AnalyticalColliderManager>();
        }
    }

    void CreateSphere()
    {
        _sphereObj = GameObject.CreatePrimitive(PrimitiveType.Sphere);
        _sphereObj.name = "TestCollider_Sphere";
        _sphereObj.transform.position = sphereStartPos;
        _sphereObj.transform.localScale = Vector3.one * sphereRadius * 2f;

        // 移除默认的Collider（我们用解析碰撞）
        // 但保留一个用于AnalyticalColliderSource自动检测
        var existingCollider = _sphereObj.GetComponent<SphereCollider>();
        if (existingCollider != null)
        {
            existingCollider.isTrigger = true; // 设为trigger避免物理干扰
        }

        // 添加解析碰撞源
        var source = _sphereObj.AddComponent<AnalyticalColliderSource>();
        source.autoDetect = true;
        source.colliderType = AnalyticalColliderType.Sphere;

        if (sphereMaterial != null)
        {
            _sphereObj.GetComponent<MeshRenderer>().material = sphereMaterial;
        }
    }

    void CreateBox()
    {
        _boxObj = GameObject.CreatePrimitive(PrimitiveType.Cube);
        _boxObj.name = "TestCollider_Box";
        _boxObj.transform.position = boxStartPos;
        _boxObj.transform.localScale = boxSize;

        // 保留BoxCollider用于自动检测
        var existingCollider = _boxObj.GetComponent<BoxCollider>();
        if (existingCollider != null)
        {
            existingCollider.isTrigger = true;
        }

        // 添加解析碰撞源
        var source = _boxObj.AddComponent<AnalyticalColliderSource>();
        source.autoDetect = true;
        source.colliderType = AnalyticalColliderType.Box;

        if (boxMaterial != null)
        {
            _boxObj.GetComponent<MeshRenderer>().material = boxMaterial;
        }
    }

    void Update()
    {
        // UpdateSphereMotion();
        // UpdateBoxMotion();
    }

    void UpdateSphereMotion()
    {
        if (_sphereObj == null) return;

        if (_spherePaused)
        {
            _spherePauseTimer += Time.deltaTime;
            if (_spherePauseTimer >= pauseDuration)
            {
                _spherePaused = false;
                _spherePauseTimer = 0f;
                _sphereMovingForward = !_sphereMovingForward;
            }
            return;
        }

        _sphereTimer += Time.deltaTime * sphereMoveSpeed;

        float t = _sphereMovingForward ? _sphereTimer : (1f - _sphereTimer);
        t = Mathf.Clamp01(t);

        // 使用SmoothStep让运动更自然
        float smoothT = t * t * (3f - 2f * t);
        _sphereObj.transform.position = Vector3.Lerp(sphereStartPos, sphereEndPos, smoothT);

        // 到达端点后停顿
        if (_sphereTimer >= 1f)
        {
            _sphereTimer = 0f;
            _spherePaused = true;
            _spherePauseTimer = 0f;
        }
    }

    void UpdateBoxMotion()
    {
        if (_boxObj == null) return;

        // 缓慢旋转
        _boxObj.transform.Rotate(Vector3.up, 30f * Time.deltaTime);

        if (_boxPaused)
        {
            _boxPauseTimer += Time.deltaTime;
            if (_boxPauseTimer >= pauseDuration)
            {
                _boxPaused = false;
                _boxPauseTimer = 0f;
                _boxMovingForward = !_boxMovingForward;
            }
            return;
        }

        _boxTimer += Time.deltaTime * boxMoveSpeed;

        float t = _boxMovingForward ? _boxTimer : (1f - _boxTimer);
        t = Mathf.Clamp01(t);

        float smoothT = t * t * (3f - 2f * t);
        _boxObj.transform.position = Vector3.Lerp(boxStartPos, boxEndPos, smoothT);

        if (_boxTimer >= 1f)
        {
            _boxTimer = 0f;
            _boxPaused = true;
            _boxPauseTimer = 0f;
        }
    }

    void OnDestroy()
    {
        if (_sphereObj != null) Destroy(_sphereObj);
        if (_boxObj != null) Destroy(_boxObj);
    }
}
