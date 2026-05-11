# DOTS-XPBD Demo 录屏脚本说明

本目录是为**毕业论文盲审 Demo 视频**准备的全套自动化录屏脚本。
所有代码均与论文实现解耦：不修改任何 `Assets/Scripts/XPBDRelated` 下的求解器/系统，只通过组合调用 `ClothRuntimeSpawner`、`SoftBodyRuntimeSpawner`、`AnalyticalColliderSource` 来搭建演示场景。

## 一句话使用流程

1. **新建场景** `Assets/Scenes/DemoRecord.unity`（File → New Scene → Empty (HDRP)/(URP) Empty 都可以，因为本工程是 URP）。
2. 场景里**新建一个空 GameObject**，命名 `DemoDirector`，挂上脚本 `DemoDirector.cs`。
3. 直接 **Play** 即可——所有灯光/相机/地板/HUD 都会程序化生成，然后自动按时间表演完 7 段，共 70 秒。
4. 用 **Unity Recorder** 录屏（见下文）。

> ⚠️ 这个场景是 **空的**，唯一一个 GameObject 就是 `DemoDirector`。**不要往里再放任何东西**——会和导演里程序化生成的内容冲突。

## 时间表（共 70 秒）

| 段 | 时长 | 内容 | HUD 标题 |
|---|---|---|---|
| 1 | 10s | 布料下落覆盖球 | 布料下落覆盖球体 |
| 2 | 10s | 布料两角悬挂自由摆动 | 布料悬挂自由摆动 |
| 3 | 10s | 布料受周期性风力（旗帜） | 布料受风飘动（旗帜） |
| 4 | 10s | 软体球下落压缩回弹 | 软体下落压缩与回弹 |
| 5 | 10s | 多软体（球/立方体）落地碰撞 | 软体多体碰撞与形变 |
| 6 | 10s | 布料盖在软体球上跨体碰撞 | 布料覆盖软体（跨体碰撞） |
| 7 | 10s | 256×256 大布料 + 3 软体球依次落下 | 大规模布料 × 多软体（性能演示） |

## 用 Unity Recorder 录视频（5.1.6 已安装）

1. 顶部菜单 `Window → General → Recorder → Recorder Window`。
2. 点 `+ Add Recorder` → `Movie`。
3. 关键设置：
   - **Source**: `Game View`
   - **Output Resolution**: `FHD - 1080p` 或 `Custom 2560×1440`（盲审视频建议 1080p 即可）
   - **Frame Rate**: `Constant 60`
   - **Recording Mode**: `Time Interval`，`Start: 0`，`End: 70`（正好 70 秒，匹配 7×10）
   - **Encoder**: `H.264 MP4` `High`
   - **Output File**: 自己选个路径，比如 `Assets/../Recordings/demo_<Take>.mp4`
4. 在 Recorder Window 点 **START RECORDING**，它会自动按 Play → 录满 70s → 停。
5. 录完去 Output File 路径取 `.mp4` 文件。

> 提示：第一次录之前，把 Game 视图的 Aspect 设成 `16:9` 或固定 `1920×1080`，避免 HUD 排版被异常长宽比拉花。

## 文件结构

```
Assets/Scripts/Demo/
├── DemoDirector.cs                # 总导演，建场景、轮播段
├── DemoSegmentBase.cs             # 段抽象基类
├── DemoCameraPreset.cs            # 相机预设（极坐标 + 起止值）
├── DemoCameraRig.cs               # 极简相机控制器（环绕/拉远/SmoothStep 插值）
├── Hud/
│   └── DemoHud.cs                 # IMGUI HUD：标题/副标题/粒子数/FPS/进度条
├── Wind/
│   └── WindForceDemo.cs           # 录屏专用伪风力（周期性脉冲叠加 ParticleVelocity）
└── Segments/
    ├── Seg01_ClothDropOnSphere.cs
    ├── Seg02_ClothHangSwing.cs
    ├── Seg03_ClothWind.cs
    ├── Seg04_SoftBodyCompress.cs
    ├── Seg05_SoftBodyCollision.cs
    ├── Seg06_CrossBodyClothOnSoft.cs
    └── Seg07_LargeScaleStress.cs
```

## 我能改什么

- **想去掉某段**：在 `DemoDirector.BuildSegments()` 里把对应那行 `_segments.Add(...)` 注释掉。
- **想换段顺序**：调换 `_segments.Add(...)` 调用顺序即可。
- **想改某段时长**：在该段类里覆写 `Duration`。
- **想改镜头**：改该段类里的 `Camera` 属性（`DemoCameraPreset.Static / Orbit / OrbitAndZoom`）。
- **想换布料/软体颜色**：在 `DemoDirector` 的 Inspector 上拖你自己的 Material 进 `ClothMaterial`/`SoftBodyMaterial`/`BallMaterial`；不指定就用默认色。
- **想关 HUD**：把 `Demo_HUD` 物体 `SetActive(false)`，或在 `DemoDirector.BuildEnvironment` 里把 HUD 那段注释掉。

## 已知约束

- 代码里没有用 Cinemachine 包，因为只需要"按段切换 + 段内插值"，自己写 5 行就够了；不影响视觉效果。
- "风力"是录屏专用伪实现（`WindForceDemo`），不写回求解器、不影响论文里的 XPBD 实现。
- 第 7 段的 `useGraphColoring=true` + `enableSelfCollision=false` 是按你提供的 256×256 = 2.77ms 性能数据反推的"稳跑 60fps"配置。
  如果你的机器不行（FPS 掉到 < 30），把 `segments / subdivision` 都降到 `127` 即可（≈16k 粒子）。
