# ext_observer · 源码解析

> 机械臂外部力矩观测器库 · 配套阅读文档
> 左边开源码，右边开这个。文件链接可点击跳转。
> 公式用 LaTeX 书写，VSCode 内置 Markdown 预览（`Ctrl+Shift+V`）原生支持渲染。

---

## 目录

1. [项目定位](#1-项目定位)
2. [核心数学：一个方程撑起整个库](#2-核心数学一个方程撑起整个库)
3. [代码地图](#3-代码地图)
4. [架构：两条继承链](#4-架构两条继承链)
5. [基础设施逐文件详解](#5-基础设施逐文件详解)
6. [七个观测器逐个详解](#6-七个观测器逐个详解)
7. [观测器选型对比](#7-观测器选型对比)
8. [测试与 MATLAB 封装](#8-测试与-matlab-封装)
9. [建议阅读顺序](#9-建议阅读顺序)
10. [阅读时的坑](#10-阅读时的坑实测确认)

---

## 1. 项目定位

估计机械臂受到的**外部力矩** $\tau_{ext}$ —— 不用任何力矩传感器，只靠关节角度 $q$、速度 $\dot q$、电机力矩 $\tau$ 和一个动力学模型。

用途是**无传感器碰撞检测**：机器人撞到人或物体时 $\tau_{ext}$ 突然跳变，超过阈值就触发停机。

与 [dynamic_calibration](../dynamic_calibration/) 是**同一篇论文的两半**：

```
dynamic_calibration  →  辨识出准确的动力学参数 (M, C, G, 摩擦, 驱动增益)
                              ↓  喂给
ext_observer         →  用这些参数实时估计 τ_ext  →  碰撞检测
```

论文：Mamedov & Mikhel, *Practical Aspects of Model-Based Collision Detection*, Frontiers in Robotics and AI, 2020。

库里实现了 **7 种观测器**，每种都有独立的理论来源（De Luca、Mohammadi、Garofalo、Hu 等）。这个库的价值在于**把这些方法放在同一套接口下做横向对比** —— 这正是论文的主题。

---

## 2. 核心数学：一个方程撑起整个库

**先读这一节。** 七个观测器全部从同一个方程出发；理解了它，剩下的只是「用什么手段解同一个问题」。

### 2.1 从动力学到动量方程

带外力的机械臂动力学：

$$
M(q)\,\ddot q + C(q,\dot q)\,\dot q + G(q) = \tau + \tau_{ext}
$$

其中 $M \in \mathbb{R}^{n\times n}$ 是惯性矩阵，$C$ 是科氏/离心矩阵，$G$ 是重力项，$\tau$ 是电机力矩，$\tau_{ext}$ 是待估计的外部力矩。

定义**广义动量**：

$$
p \;=\; M(q)\,\dot q
$$

对时间求导（注意 $M$ 依赖 $q$，$q$ 依赖 $t$，所以是乘积求导）：

$$
\dot p \;=\; \dot M\,\dot q + M\,\ddot q
$$

从动力学方程解出 $M\ddot q = \tau + \tau_{ext} - C\dot q - G$，代入：

$$
\dot p \;=\; \dot M\,\dot q + \tau + \tau_{ext} - C\dot q - G
$$

### 2.2 关键一步：斜对称性

这里用上机械臂动力学最重要的结构性质。矩阵 $\dot M - 2C$ 是**反对称**的，即 $x^T(\dot M - 2C)x = 0\ \forall x$。由此可推出（在 Christoffel 符号定义的 $C$ 下）：

$$
\boxed{\;\dot M \;=\; C + C^T\;}
$$

代入上式：

$$
\dot p = (C + C^T)\dot q + \tau + \tau_{ext} - C\dot q - G
       = C^T\dot q - G + \tau + \tau_{ext}
$$

定义 $\beta$（代码里就叫 `beta`）：

$$
\beta(q,\dot q) \;\triangleq\; G(q) - C(q,\dot q)^T\,\dot q
$$

得到**全库的核心方程**：

$$
\boxed{\;\dot p \;=\; \tau - \beta + \tau_{ext}\;}
$$

> 这个方程会在后面每个观测器的推导里反复出现。建议现在就在纸上推一遍。

### 2.3 为什么不能直接求解

移项就得到答案：

$$
\tau_{ext} = \dot p - \tau + \beta
$$

**但这行不通。** $\dot p$ 需要数值微分：$p = M(q)\dot q$ 里的 $\dot q$ 本身就带噪（往往还是编码器差分出来的），再微分一次，噪声被放大到完全淹没信号。频域上看，微分器的幅频响应是 $|j\omega| = \omega$ —— 频率越高增益越大，而噪声恰恰集中在高频。

而且 $\dot p$ 还需要 $\dot M$，解析算很麻烦。

**所以：七个观测器 = 七种规避微分 $p$ 的手段。** 这是理解整个库的钥匙：


| 观测器      | 规避微分的手段                          | 代价                   |
| ----------- | --------------------------------------- | ---------------------- |
| Momentum    | 积分反馈，$r$ 一阶低通逼近 $\tau_{ext}$ | 需整定增益$K_O$        |
| Disturbance | $M^{-1}$ 加权的一阶滤波 + 隐式离散      | 每步 2 次矩阵求逆      |
| SlidingMode | super-twisting 二阶滑模，有限时间收敛   | 4 组增益，可能抖振     |
| DKalman     | 把$\tau_{ext}$ 建成状态变量，卡尔曼估计 | 需整定$Q,R$            |
| DKalmanExp  | 同上 + 精确离散化 + 位形相关的$R$       | 最贵                   |
| FDyn        | 滤波器代数恒等式把微分**吸收**掉        | 无反馈，不抑制模型误差 |
| FRange      | FDyn + 模型不确定性带                   | 同上                   |

### 2.4 记号约定


| 数学符号                | 代码变量  | 含义                                     |
| ----------------------- | --------- | ---------------------------------------- |
| $q,\ \dot q$            | `q`, `qd` | 关节位置、速度                           |
| $\ddot q$               | `q2d`     | 关节加速度（观测器**不需要**它）         |
| $\tau$                  | `tau`     | 测量到的电机力矩（**已换算**，不是电流） |
| $p = M\dot q$           | `p`       | 广义动量                                 |
| $\beta = G - C^T\dot q$ | `beta`    | 核心方程里的组合项                       |
| $\tau_{ext}$            | 返回值    | 待估计的外部力矩                         |
| $\Delta t$              | `dt`      | 距上次调用的时间步                       |
| $n$                     | `jointNo` | 关节数                                   |

> ⚠️ **`tau` 必须是力矩，不是电流。** 电流 → 力矩的驱动增益 $K$ 由调用方负责（`test.m` 里的 `K.*cur(i,:)`），而那个 $K$ 正是 `dynamic_calibration` 的 `estimate_drive_gains` 辨识出来的。**两个项目在这里对接。**

---

## 3. 代码地图

整个库只有 **28 个文件**，**零第三方依赖**（除 Eigen）。可以全部读完。

```
ext_observer/
├── lib/                          ← 库本体，头文件即实现
│   ├── external_observer.h       ← 抽象基类（先读这个）
│   ├── robot_dynamics_rnea.cpp   ← RNEA 数值技巧（★全库最精彩）
│   │
│   ├── momentum_observer.h            ┐
│   ├── disturbance_observer.h         │
│   ├── sliding_mode_observer.h        ├ 7 个观测器（M/C/G 版）
│   ├── disturbance_kalman_filter.h    │
│   ├── disturbance_kalman_filter_exp.h│
│   ├── filtered_dyn_observer.h        │
│   ├── filtered_range_observer.h      ┘
│   │
│   ├── momentum_observer_rnea.h        ┐
│   ├── disturbance_observer_rnea.h     │
│   ├── sliding_mode_observer_rnea.h    ├ 5 个观测器（RNEA 版）
│   ├── disturbance_kalman_filter_rnea.h│
│   ├── filtered_dyn_observer_rnea.h    ┘
│   │
│   ├── iir_filter.h              ← 数字滤波器（双线性变换）
│   ├── kalman_filter.h           ← 离散卡尔曼
│   ├── kalman_filter_continous.h ┐
│   └── kalman_filter_continous.cpp┘ 连续系统卡尔曼（Van Loan 离散化）
│
├── tests/                        ← 示例 + 人眼验证（非自动化测试）
│   ├── double_link.h             ← 2 连杆玩具模型（M/C/G 版）
│   ├── double_link_rnea.h        ← 同上（RNEA 版）
│   ├── test.cpp                  ← M/C/G 观测器示例
│   ├── test_rnea.cpp             ← RNEA 观测器示例
│   ├── test_kalman.cpp           ← 卡尔曼滤波器单独验证
│   ├── observers.cpp/.h          ← ⚠️ matlab/observers.cpp 的过时副本
│   ├── Makefile
│   └── *.gnuplot                 ← 画图脚本
│
└── matlab/                       ← MATLAB 封装（C 接口 + dlopen）
    ├── observers.cpp/.h          ← extern "C" 封装层
    ├── test.m                    ← 七种观测器横向对比（★论文实验）
    ├── test_lib.cpp              ← 从 C++ 调用 .so 的示例
    └── Makefile
```

### 3.1 header-only 的结构性后果

几乎所有代码都在 `.h` 里**包括函数体**，且函数**既不是 `inline` 也不是模板**。只有 3 个 `.cpp`：`robot_dynamics_rnea.cpp`、`kalman_filter_continous.cpp` 和封装层。

这带来一个硬约束：

> **同一个观测器头文件不能被两个编译单元同时 include**，否则违反 ODR（One Definition Rule），链接期重复定义报错。

所以 [tests/test.cpp](tests/test.cpp) 用 `#ifdef` 一次只启用一个观测器 —— **这不是为了方便切换，是结构性约束**。而 [matlab/observers.cpp](matlab/observers.cpp) 能一次 include 全部 7 个，是因为它是整个 `.so` 里唯一的编译单元。

如果你要在自己的项目里用多个观测器且分多个 `.cpp`，必须先把这些头文件改成 `inline` 或拆出 `.cpp`。

---

## 4. 架构：两条继承链

这是读这个库最先要建立的心智模型。**每个观测器都有两个版本**，源头是「机器人动力学怎么提供」这个问题有两种答案。

```
        RobotDynamicsBase                          ExternalObserverBase
        (getFriction, jointNo)                (getExternalTorque, reset, type)
           ↙            ↘                              ↙              ↘
  RobotDynamics    RobotDynamicsRnea         ExternalObserver    ExternalObserverRnea
  (getM/getC/getG) (rnea + getM/tranCqd)    (持有RobotDynamics*)(持有RobotDynamicsRnea*)
        ↑                  ↑                          ↑                   ↑
    DoubleLink      DoubleLinkRnea            MomentumObserver   MomentumObserverRnea
   （你的机器人）    （你的机器人）              DisturbanceObserver       ...
                                                    ...
```

**左链 `RobotDynamics`**：你能显式给出 $M(q)$、$C(q,\dot q)$、$G(q)$。适合自己推导过动力学的场合 —— 比如 `dynamic_calibration` 生成的 `M_mtrx_fcn.m` / `C_mtrx_fcn.m` / `G_vctr_fcn.m`。

**右链 `RobotDynamicsRnea`**：你只有一个 RNEA 黑盒 $f(q,\dot q,\ddot q) \to \tau$。这是绝大多数动力学库（Orocos KDL、Pinocchio、RBDL）提供的形式。观测器需要的 $M$、$C^T\dot q$ **全部靠巧妙地调用这个黑盒凑出来** —— 见 §5.2。

### 4.1 使用方式

实现一个接口类，把指针交给观测器：

```cpp
class MyRobot : public RobotDynamics {
  MatrixJ getM(VectorJ& q) override               { /* ... */ }
  MatrixJ getC(VectorJ& q, VectorJ& qd) override  { /* ... */ }
  VectorJ getG(VectorJ& q) override               { /* ... */ }
  VectorJ getFriction(VectorJ& qd) override       { /* ... */ }
  int jointNo() override { return 6; }
};

MyRobot robot;
VectorJ k(6);  k << 30,30,30,30,30,30;
MomentumObserver obs(&robot, k);

// 控制循环里
VectorJ ext = obs.getExternalTorque(q, qd, tau, dt);
```

**统一入口只有一个方法**（[lib/external_observer.h:162](lib/external_observer.h#L162)）：

```cpp
virtual VectorJ getExternalTorque(VectorJ& q, VectorJ& qd, VectorJ& tau, double dt) = 0;
```

七个观测器可以随意互换 —— 这就是论文能做横向对比的基础。

### 4.2 观测器 ID

每个观测器有个整型 ID（`objType`），用于运行时类型识别：


| ID | 观测器                | RNEA 版 ID |
| -- | --------------------- | ---------- |
| 1  | `DKalmanObserver`     | 21         |
| 2  | `FDynObserver`        | 22         |
| 3  | `SlidingModeObserver` | 23         |
| 4  | `DisturbanceObserver` | 24         |
| 5  | `MomentumObserver`    | 25         |
| 6  | `FRangeObserver`      | —         |
| 7  | `DKalmanObserverExp`  | —         |

**规律：RNEA 版 $=$ 基础版 $+\ 20$。** `FRangeObserver` 和 `DKalmanObserverExp` 没有 RNEA 版本。

这个 ID 存在的唯一原因是 [matlab/observers.cpp](matlab/observers.cpp) 的 `freeAll()` 要在只有基类指针时 `delete` 出正确的派生类型（§8.2 会说明那其实是个设计缺陷）。

### 4.3 通用状态机：`isRun`

所有观测器都是有状态的（要么积分、要么滤波、要么 KF）。共用同一个首次调用模式：

```cpp
if(isRun) {
    // 正常更新
} else {
    // 用当前测量初始化内部状态
    isRun = true;
}
```

**所有观测器首次调用一律返回 $0$** —— 我逐个追过初始化路径，这个不变量是成立的（§6 每个小节都有「首次调用路径」的追踪）。

`reset()` 把 `isRun` 置 false，下次调用重新初始化：

```cpp
inline void reset() noexcept { isRun = false; }
```

> **机器人重新上电、或长时间暂停后必须调用 `reset()`**，否则积分器/滤波器里的陈旧状态会产生假报警。

---

## 5. 基础设施逐文件详解

### 5.1 [lib/external_observer.h](lib/external_observer.h) —— 类型与抽象接口

#### 类型定义

```cpp
const double GRAVITY = 9.81;
typedef Eigen::MatrixXd MatrixJ;   // 动态大小矩阵
typedef Eigen::VectorXd VectorJ;   // 动态大小向量
```

`MatrixXd` / `VectorXd` 是**动态大小**的 Eigen 类型 —— 意味着堆分配。对硬实时控制这不是最优：固定大小的 `Matrix<double,6,6>` 会栈分配、能向量化、无分配开销。这里换来的是「一套代码支持任意关节数」。**这个权衡贯穿全库**。

> **历史包袱**：这两个 typedef 是 2024-12 的 `675b5da Update code style` 提交引入的，之前叫 `Vector` / `Matrix`。这次改名留下了几处未更新的死代码，见 §10.2。

#### `RobotDynamicsBase` —— 两条链的公共部分

```cpp
class RobotDynamicsBase {
public:
  virtual ~RobotDynamicsBase() = default;
  virtual VectorJ getFriction(VectorJ& qd) = 0;   // 摩擦模型
  virtual int jointNo() = 0;                      // 关节数
};
```

只有两个纯虚函数。**摩擦对两条链都是必须的** —— 因为观测器要从测量力矩里扣掉摩擦，否则摩擦会被误判成外力。

`getFriction` 只依赖 $\dot q$，说明库假设的是**静态摩擦模型**（如 $\tau_f = F_v\dot q + F_c\,\mathrm{sign}(\dot q) + F_0$，正是 `dynamic_calibration` 的 `frictionRegressor` 辨识的那个），不支持 LuGre 这类有内部状态的动态摩擦。

#### `RobotDynamics` —— 显式矩阵链

```cpp
class RobotDynamics : public RobotDynamicsBase {
public:
  virtual MatrixJ getM(VectorJ& q) = 0;
  virtual MatrixJ getC(VectorJ& q, VectorJ& qd) = 0;
  virtual VectorJ getG(VectorJ& q) = 0;
};
```

三个纯虚函数，语义直白。注意 **`getC` 返回的是矩阵 $C$ 本身**（不是乘积 $C\dot q$）—— 因为观测器需要 $C^T\dot q$，必须拿到矩阵才能转置。这是左链相对右链的关键优势。

#### `RobotDynamicsRnea` —— RNEA 链

```cpp
class RobotDynamicsRnea : public RobotDynamicsBase {
public:
  RobotDynamicsRnea();
  virtual VectorJ rnea(VectorJ& q, VectorJ& qd, VectorJ& q2d, double g = 0) = 0;
  inline void setDelta(double d) noexcept { delta = d; }
  VectorJ tranCqd(VectorJ& q, VectorJ& qd);   // ← 已实现，非纯虚
  MatrixJ getM(VectorJ& q);                   // ← 已实现，非纯虚

protected:
  VectorJ _qext, _p0, _zero, _sum;   /**< Temporary variables, avoid memory reallocation. */
  MatrixJ _M;
  double delta = 1E-7;               /**< 数值微分步长 */
};
```

**只有 `rnea` 是纯虚的** —— 你只需实现这一个，`getM` 和 `tranCqd` 由基类用 `rnea` 凑出来（§5.2）。这是很聪明的设计：把「RNEA → 观测器需要的量」这个转换的复杂度一次性封装在基类里。

注意 `double g = 0` —— **重力默认关闭**。读 `*_rnea.h` 时必须盯紧哪次调用传了 `GRAVITY`。

#### 临时变量成员：一个打了折扣的优化

`_qext, _p0, _zero, _sum, _M` 都是成员而非局部变量，注释写着 `avoid memory reallocation`。这是实时控制代码的标准模式 —— 避免每个控制周期堆分配。全库的观测器都这么做（`sum`、`r`、`p`、`beta`、`torque`、`tprev`…）。

> **但这个努力被接口设计部分抵消了**：`getM()` 返回 `MatrixJ`（**按值**），于是 `dyn->getM(q) * qd` 这一行里，`getM` 的返回值是临时对象，乘法结果又是一个临时对象。真要做到零分配，接口得改成 `getM(VectorJ& q, MatrixJ& out)`。
>
> 读代码时不必被 `// Temporary objects` 注释误导 —— 它们**减少**了分配，但没有**消除**。

#### `ExternalObserverBase`

```cpp
class ExternalObserverBase {
public:
  ExternalObserverBase() : jointNo(0), objType(-1), isRun(false) {}
  virtual ~ExternalObserverBase() = default;      // ★ 虚析构（§10.8 会用到）
  virtual VectorJ getExternalTorque(VectorJ& q, VectorJ& qd, VectorJ& tau, double dt) = 0;
  inline void reset() noexcept { isRun = false; }
  inline int type() const noexcept { return objType; }
protected:
  int jointNo;    /**< 关节数 */
  int objType;    /**< 观测器类型 ID */
  bool isRun;     /**< 首次调用标志 */
};
```

**`virtual ~ExternalObserverBase() = default;` 很重要** —— 有了它，`delete` 一个基类指针就能正确调用派生类析构。§10.8 会说明 `matlab/observers.cpp` 里那个巨大的 `switch` 因此完全多余。

#### `ExternalObserver` / `ExternalObserverRnea`

```cpp
class ExternalObserver : public ExternalObserverBase {
public:
  ExternalObserver(RobotDynamics *rd, int type)
  : ExternalObserverBase(), dyn(rd)
  {
    jointNo = rd->jointNo();    // ← 构造时缓存关节数
    objType = type;
  }
protected:
  RobotDynamics *dyn;
};
```

两个类除了持有的指针类型不同，**完全一样**。构造时把 `jointNo` 从机器人对象里拷出来缓存 —— 这样所有派生类的初始化列表里就能直接用 `VectorJ(jointNo)` 分配临时量。

注意**基类成员 `jointNo` 是在构造函数体里赋值的**，而派生类的初始化列表（如 `sum(VectorJ(jointNo))`）在基类构造**完成之后**才执行 —— 所以顺序是安全的。这是个容易忽略但必须成立的细节。

> `dyn` 是**裸指针，不拥有所有权**。观测器不会 delete 它。调用方必须保证机器人对象的生命周期长于观测器 —— `matlab/observers.cpp` 用 `static DoubleLink robot;` 全局对象来保证这点。

---

### 5.2 ★ [lib/robot_dynamics_rnea.cpp](lib/robot_dynamics_rnea.cpp) —— 全库最精彩的 55 行

**问题**：观测器需要 $M(q)$ 和 $C^T\dot q$，但 RNEA 只给你一个黑盒：

$$
\texttt{rnea}(q, \dot q, \ddot q, g) \;=\; M(q)\,\ddot q + C(q,\dot q)\,\dot q + G(q, g)
$$

#### RNEA 调用约定速查表

**这张表是读所有 `*_rnea.h` 文件的钥匙**，见到 `rnea(...)` 就查它：


| 调用                      | 得到                  | 原理                                                                   |
| ------------------------- | --------------------- | ---------------------------------------------------------------------- |
| `rnea(q, 0, qd)`          | $M(q)\,\dot q$        | 把$\dot q$ 塞进**加速度**位置！$\dot q{=}0$ 消掉 $C$，$g{=}0$ 消掉 $G$ |
| `rnea(q, 0, 0, GRAVITY)`  | $G(q)$                | $\ddot q{=}0,\ \dot q{=}0$，只剩重力                                   |
| `rnea(q, qd, 0)`          | $C(q,\dot q)\,\dot q$ | $\ddot q{=}0$ 消掉 $M$，$g{=}0$ 消掉 $G$                               |
| `rnea(q, qd, 0, GRAVITY)` | $C\dot q + G$         | 上面两个合并                                                           |
| `rnea(q, 0, e_i)`         | $M$ 的第 $i$ 列       | $e_i$ 是单位向量                                                       |

第一行是最容易看晕的：

```cpp
p = dyn->rnea(q, zero, qd);     // M * qd
```

**速度向量被当作加速度传进去。** 因为 RNEA 对 $\ddot q$ 是线性的：

$$
\texttt{rnea}(q, 0, x, 0) = M(q)\,x \qquad \forall x
$$

取 $x = \dot q$ 就白嫖到了动量。这个技巧在全库出现十几次。

#### `getM()` —— 逐列构造惯性矩阵

```cpp
MatrixJ RobotDynamicsRnea::getM(VectorJ& q)
{
  int N = jointNo();
  _zero.resize(N);  _zero.setZero();
  _M.resize(N, N);
  _qext.resize(N);

  for(int i = 0; i < N; i++) {
    _qext.setZero();                  // 复用 _qext 当作单位向量 e_i
    _qext(i) = 1;
    _p0 = rnea(q, _zero, _qext);      // 复用 _p0 存放矩阵的一列
    for(int j = 0; j < N; j++)
      _M(j, i) = _p0(j);
  }
  return _M;
}
```

原理：

$$
\texttt{rnea}(q, 0, e_i, 0) = M(q)\,e_i = M_{:,i}
$$

即 $M$ 的第 $i$ 列。跑 $N$ 次就拼出整个矩阵。

**逐行注意点**：

- `_qext` 和 `_p0` 都被**复用**了（注释明说 `use _qext to choose column` / `use _p0 to save matrix column`），它们原本的语义（"扩展的 q"、"初始动量"）在这里完全不成立。这是全库变量复用的典型样本。
- 内层 `for(j)` 手动逐元素拷贝，其实 `_M.col(i) = _p0;` 一行就够 —— Eigen 会做得更快。
- **成本：$N$ 次 RNEA 调用**（6 轴 = 6 次）。

#### `tranCqd()` —— 数值求 $C^T\dot q$（★核心技巧）

$C^T\dot q$ 没法直接从 RNEA 拿到 —— RNEA 只给 $C\dot q$，**转置不出来**。这里用斜对称恒等式（§2.2）绕开：

$$
\dot M = C + C^T \quad\Longrightarrow\quad C^T\dot q = \dot M\,\dot q - C\,\dot q
$$

$C\dot q$ 好办（查表第三行）。问题变成怎么求 $\dot M\dot q$。用链式法则：

$$
\dot M = \sum_{i=1}^{n} \frac{\partial M}{\partial q_i}\,\dot q_i
\quad\Longrightarrow\quad
\dot M\,\dot q = \sum_{i=1}^{n} \left(\frac{\partial M}{\partial q_i}\,\dot q\right)\dot q_i
$$

而括号里的 $\dfrac{\partial M}{\partial q_i}\dot q$ 可以用**前向差分**逼近，且每一项恰好是一次 RNEA 调用：

$$
\frac{\partial M}{\partial q_i}\,\dot q \;\approx\;
\frac{M(q + \delta e_i)\,\dot q - M(q)\,\dot q}{\delta}
= \frac{\texttt{rnea}(q + \delta e_i,\, 0,\, \dot q) - \texttt{rnea}(q,\, 0,\, \dot q)}{\delta}
$$

**逐行对照**：

```cpp
VectorJ RobotDynamicsRnea::tranCqd(VectorJ& q, VectorJ& qd)
{
  int N = jointNo();
  _zero.resize(N);  _zero.setZero();
  _p0 = rnea(q, _zero, qd);          // ① _p0 = M(q)·q̇     ← 差分的基准点
  _sum.resize(N);   _sum.setZero();

  for(int i = 0; i < N; i++) {
    _qext = q;                       // ② 拷贝 q
    _qext(i) += delta;               //    q + δe_i
    _sum += (rnea(_qext, _zero, qd) - _p0) * (qd(i) / delta);
    //       └── M(q+δe_i)·q̇ ──┘   └M·q̇┘    └── ×q̇_i/δ ──┘
  }
  // ③ 此时 _sum = Ṁ·q̇

  _sum -= rnea(q, qd, _zero);        // ④ 减去 C·q̇
  // ⑤ 此时 _sum = Ṁ·q̇ - C·q̇ = Cᵀ·q̇

  return _sum;
}
```


| 步骤 | 数学                                                                              | 说明                                            |
| ---- | --------------------------------------------------------------------------------- | ----------------------------------------------- |
| ①   | $p_0 = M(q)\dot q$                                                                | 差分基准，循环里$N$ 次都用它，只算一次          |
| ②   | $q + \delta e_i$                                                                  | 只扰动第$i$ 个关节                              |
| ③   | $\sum_i \frac{M(q+\delta e_i)\dot q - M(q)\dot q}{\delta}\dot q_i = \dot M\dot q$ | 除法和乘$\dot q_i$ 合并成一次 `* (qd(i)/delta)` |
| ④   | $-\,C\dot q$                                                                      | `rnea(q, qd, 0)` 查表第三行                     |
| ⑤   | $C^T\dot q$                                                                       | ✓                                              |

> ⚠️ **注释 `// M'*qd - C*qd` 里的 `M'` 指的是 $\dot M = dM/dt$，不是转置 $M^T$。** 在一个到处是 `.transpose()` 的文件里写 `M'`，极易误读。这是我读这个库时卡住最久的一行。

**成本：$N + 2$ 次 RNEA 调用**（1 次基准 + $N$ 次循环 + 1 次 $C\dot q$；6 轴 = 8 次）。

#### 差分步长 `delta = 1E-7`

前向差分的总误差 $\approx \underbrace{O(\delta)}_{\text{截断}} + \underbrace{O(\epsilon/\delta)}_{\text{舍入}}$，其中 $\epsilon \approx 2.2\times10^{-16}$ 是双精度机器精度。最优步长：

$$
\delta^* \approx \sqrt{\epsilon} \approx 1.5\times 10^{-8}
$$

代码取 $10^{-7}$，同一量级，合理。可以用 `setDelta()` 在运行时调整。

> 若改用**中心差分** $\frac{M(q+\delta e_i)\dot q - M(q-\delta e_i)\dot q}{2\delta}$，精度提升到 $O(\delta^2)$，但 RNEA 调用次数翻倍到 $2N+1$。作者选了前向差分 —— 实时性优先。

#### RNEA 版本的性能代价

这是选择实现路线时的核心权衡：


| 观测器                    | 每步 RNEA 调用次数                                                            |
| ------------------------- | ----------------------------------------------------------------------------- |
| `MomentumObserverRnea`    | $1_{(M\dot q)} + 1_{(G)} + 8_{(\texttt{tranCqd})} = \mathbf{10}$              |
| `FDynObserverRnea`        | $1 + 1 + 8 = \mathbf{10}$                                                     |
| `SlidingModeObserverRnea` | $1 + 1 + 8 = \mathbf{10}$                                                     |
| `DKalmanObserverRnea`     | $1 + 1 + 8 = \mathbf{10}$                                                     |
| `DisturbanceObserverRnea` | $6_{(\texttt{getM})} + 1_{(C\dot q + G)} = \mathbf{7}$ + 两次 $6\times6$ 求逆 |

> **如果你有解析的 $M/C/G$（比如 `dynamic_calibration` 生成的），一定用非 RNEA 版本** —— 快约一个数量级。RNEA 版本存在的意义是兼容只提供 RNEA 的第三方动力学库。

#### `// TODO: call it once`

`tranCqd()` 和 `getM()` 开头都有这个注释，指的是：

```cpp
int N = jointNo();          // 虚函数调用
_zero.resize(N);            // 每次都 resize
_zero.setZero();            // 每次都清零
```

这些本可以在构造时做一次。是作者自己标记的优化点，不影响正确性。（`resize` 到相同尺寸时 Eigen 不会重分配，所以实际开销主要是 `setZero` 和虚函数调用。）

---

### 5.3 [lib/iir_filter.h](lib/iir_filter.h) —— 双线性变换

五个滤波器类，但**只有两个被观测器使用**：


| 类                  | 传递函数                             | 谁在用                                     |
| ------------------- | ------------------------------------ | ------------------------------------------ |
| `FilterF1`          | $H(s) = \dfrac{\omega}{s+\omega}$    | `FDynObserver`, `FRangeObserver` ✓        |
| `FilterF2`          | $H(s) = \dfrac{-\omega^2}{s+\omega}$ | 同上 ✓                                    |
| `FilterButterworth` | 二阶 Butterworth                     | ✗**死代码**                               |
| `FilterLowPass`     | 一阶低通                             | ✗**死代码**（与 `FilterF1` 系数完全相同） |
| `FilterHighPass`    | 一阶高通                             | ✗**死代码**                               |

#### 抽象基类

```cpp
class FilterIIR {
public:
  virtual Eigen::VectorXd filt(Eigen::VectorXd& x) = 0;
  virtual void update(double cutOff, double sampTime) = 0;
};
```

注意**没有虚析构函数** —— 但因为所有滤波器都是观测器的**值成员**（不是 `new` 出来的指针），never deleted polymorphically，所以实际不会出问题。

#### 预畸变双线性变换 —— 完整推导

所有滤波器都用同一个模式：

```cpp
double omega = tan(cutOff * sampTime * 0.5);   // ← 预畸变 (prewarping)
```

**双线性变换**（Tustin 法）把 $s$ 域映射到 $z$ 域：

$$
s \;\longleftarrow\; \frac{2}{T}\cdot\frac{1 - z^{-1}}{1 + z^{-1}}
$$

它的频率映射关系是非线性的：模拟频率 $\omega_a$ 与数字频率 $\omega_d$ 满足

$$
\omega_a = \frac{2}{T}\tan\frac{\omega_d T}{2}
$$

高频段被压缩，导致截止频率偏移。**预畸变**就是反向补偿：把目标 $\omega$ 先扭曲成 $\frac{2}{T}\tan\frac{\omega T}{2}$ 再代入，保证变换后截止频率精确落在期望位置。

代码里的 `omega` 变量 $= \tan\frac{\omega T}{2}$，记作 $\Omega$。

#### `FilterF1` 逐行推导

$$
H(s) = \frac{\omega}{s + \omega}
$$

代入双线性变换（用预畸变后的 $\omega \to \frac{2}{T}\Omega$）：

$$
H(z) = \frac{\frac{2}{T}\Omega}{\frac{2}{T}\frac{1-z^{-1}}{1+z^{-1}} + \frac{2}{T}\Omega}
     = \frac{\Omega(1+z^{-1})}{(1-z^{-1}) + \Omega(1+z^{-1})}
     = \frac{\Omega(1+z^{-1})}{(1+\Omega) + (\Omega-1)z^{-1}}
$$

化成差分方程：

$$
\left[(1+\Omega) + (\Omega-1)z^{-1}\right] y = \Omega(1+z^{-1})\,x
$$

$$
\boxed{\;y[k] = \underbrace{\frac{1-\Omega}{1+\Omega}}_{k_1}\,y[k{-}1] + \underbrace{\frac{\Omega}{1+\Omega}}_{k_2}\big(x[k] + x[k{-}1]\big)\;}
$$

对应实现，一一吻合：

```cpp
void FilterF1::update(double cutOff, double sampTime)
{
  double omega = tan(cutOff*sampTime*0.5);   // Ω
  k1 = (1-omega) / (1+omega);                // (1-Ω)/(1+Ω)
  k2 = omega / (1 + omega);                  // Ω/(1+Ω)
  cut = cutOff;                              // ← 保存原始 ω，getOmega() 要用
}

Eigen::VectorXd FilterF1::filt(Eigen::VectorXd& x)
{
  y1 = k1*y1 + k2*(x + x1);   // y[k] = k1·y[k-1] + k2·(x[k] + x[k-1])
  x1 = x;                     // 更新历史
  return y1;
}
```

状态变量只有两个：`x1` = $x[k{-}1]$，`y1` = $y[k{-}1]$。一阶 IIR 的最小状态。

#### `FilterF2` 逐行推导

$$
H(s) = \frac{-\omega^2}{s+\omega}
$$

同样代入，注意分子的 $\omega^2$ 也要用预畸变值 $\left(\frac{2}{T}\Omega\right)^2$。令 $f_2 = \frac{2}{T}$：

$$
H(z) = \frac{-(f_2\Omega)^2}{f_2\frac{1-z^{-1}}{1+z^{-1}} + f_2\Omega}
     = \frac{-(f_2\Omega)^2(1+z^{-1})}{f_2\left[(1-z^{-1}) + \Omega(1+z^{-1})\right]}
     = \frac{-f_2\Omega^2(1+z^{-1})}{(1+\Omega) + (\Omega-1)z^{-1}}
$$

$$
\boxed{\;y[k] = \underbrace{\frac{1-\Omega}{1+\Omega}}_{k_1}y[k{-}1] + \underbrace{\frac{-f_2\Omega^2}{1+\Omega}}_{k_2}\big(x[k]+x[k{-}1]\big)\;}
$$

```cpp
void FilterF2::update(double cutOff, double sampTime)
{
  cut = cutOff;
  f2 = 2/sampTime;                            // f₂ = 2/T
  omega = tan(cutOff*sampTime*0.5);           // Ω
  k1 = (1-omega)/(1+omega);                   // 与 F1 相同
  k2 = -f2*omega*omega/(1+omega);             // -f₂Ω²/(1+Ω)
}
```

`filt()` 的形式和 `FilterF1` **一模一样**（`y1 = k1*y1 + k2*(x + x1)`）—— 只有系数不同。两个滤波器的差异全部封装在 `k2` 里。

#### 变步长版本

```cpp
Eigen::VectorXd FilterF1::filt(Eigen::VectorXd& x, double dt)
{
  double omega = tan(cut*dt*0.5);     // 用当前 dt 重算 Ω
  k1 = (1-omega) / (1+omega);
  k2 = omega / (1 + omega);
  return filt(x);                     // 复用定步长实现
}
```

**观测器用的是这个版本** —— 因为真实控制循环的 `dt` 会抖动。这也是 `cut`（原始截止频率）必须被保存的原因：定步长 `update()` 算完系数就可以丢掉 $\omega$，变步长版本每次都要重算。

代价是每步几次 `tan()` 调用。对 1 kHz 控制环这可能不可忽略。

> `FilterButterworth` **没有**变步长版本 —— 这也从侧面印证它没被观测器用（观测器都需要变步长）。

#### `set()` —— 稳态初始化

```cpp
inline void FilterF1::set(Eigen::VectorXd& x0) { x1 = x0; y1 = x0; }
inline void FilterF2::set(Eigen::VectorXd& x0) { x1 = x0; y1 = -f2*omega*x0; }
```

把滤波器状态直接设成对应输入的**稳态输出**，消除启动瞬态。直流增益：

$$
H_{F1}(0) = \frac{\omega}{0+\omega} = 1 \quad\Longrightarrow\quad y_1 = x_0 \;\checkmark
$$

$$
H_{F2}(0) = \frac{-\omega^2}{0+\omega} = -\omega \quad\Longrightarrow\quad y_1 = -\omega\, x_0 \;\checkmark
$$

代码里 `-f2*omega*x0` 中 $f_2\Omega = \frac{2}{T}\tan\frac{\omega T}{2} \approx \omega$（小 $\omega T$ 时），正是预畸变后的 $\omega$。**这个细节做得很到位** —— 用的是离散化后的等效 $\omega$ 而不是原始 $\omega$，两者一致。

#### `getOmega()` 的命名陷阱

```cpp
inline double getOmega() const noexcept { return cut; }
```

返回的是 `cut`（**截止频率** $\omega$），**不是**内部那个 `omega` 成员（$\Omega = \tan\frac{\omega T}{2}$）。同一个类里有两个都叫 "omega" 的概念，getter 返回的是**没有** `omega` 这个名字的那个。极易误读。

但 `FDynObserver` 需要的确实是 $\omega$ —— §6.6 的推导会说明为什么。

#### ⚠️ 单位陷阱

见 §10.3：观测器构造函数的参数名 `cutOffHz` / `sampHz` **都在撒谎**，实际单位是 rad/s 和**秒**。

---

### 5.4 [lib/kalman_filter.h](lib/kalman_filter.h) —— 离散卡尔曼

#### 标准方程

系统模型：

$$
X_k = A X_{k-1} + B u_k + w_k, \qquad w_k \sim \mathcal{N}(0, Q)
$$

$$
y_k = C X_k + v_k, \qquad v_k \sim \mathcal{N}(0, R)
$$

**预测步**：

$$
\hat X_{k|k-1} = A\hat X_{k-1|k-1} + Bu_k
$$

$$
P_{k|k-1} = A P_{k-1|k-1} A^T + Q
$$

**更新步**：

$$
S_k = C P_{k|k-1} C^T + R \qquad \text{(新息协方差，代码里叫 \texttt{Y})}
$$

$$
K_k = P_{k|k-1} C^T S_k^{-1}
$$

$$
\hat X_{k|k} = \hat X_{k|k-1} + K_k\big(y_k - C\hat X_{k|k-1}\big)
$$

$$
P_{k|k} = (I - K_k C)\,P_{k|k-1}
$$

#### 逐行对照

```cpp
Eigen::VectorXd KalmanFilter::step(Eigen::VectorXd& u, Eigen::VectorXd& y)
{
  // predict
  X = A * X + B * u;                       // X̂(k|k-1)
  P = A * P * A.transpose() + Q;           // P(k|k-1)
  // update
  Y = C * (P * C.transpose()) + R;         // S = CPCᵀ + R      ← 变量名 Y 有歧义
  K = P * (C.transpose() * Y.inverse());   // K = PCᵀS⁻¹
  X += K * (y - C * X);                    // X̂(k|k)，注意此处 X 已是先验
  P = (I - K * C) * P;                     // P(k|k)
  return X;
}
```

**逐行注意点**：

- **`Y` 这个成员名有歧义** —— 它是**新息协方差** $S$，不是观测量。而函数参数里的 `y` 才是观测量。大小写区分两个完全不同的东西，是个可读性地雷。
- `X += K * (y - C * X)` —— 这里的 `C * X` 用的是**刚更新过的先验** $\hat X_{k|k-1}$，正确。
- `Y.inverse()` —— 显式求逆。对小矩阵（$n\times n$，$n\le 6$）可以接受；数值上更稳的做法是解线性方程组（`Y.ldlt().solve(...)`）。

#### 协方差更新用的是简化形式

$$
P_{k|k} = (I - K_kC)P_{k|k-1}
$$

这是**简化形式**，不是 **Joseph 形式**：

$$
P_{k|k} = (I - K_kC)P_{k|k-1}(I-K_kC)^T + K_kRK_k^T
$$

简化形式更便宜，但**只在 $K$ 恰好是最优增益时才成立**；数值误差累积会让 $P$ 逐渐失去对称性和正定性。Joseph 形式对任意 $K$ 都保持对称正定。

对短时间运行没问题，**长时间跑（数小时）可能需要换成 Joseph 形式**，或周期性地做 `P = 0.5*(P + P.transpose())` 强制对称化。

#### `reset()` 把 $P$ 置零

```cpp
void KalmanFilter::reset(Eigen::VectorXd& x0)
{
  X = x0;
  P = Eigen::MatrixXd::Zero(nx, nx);    // ← P₀ = 0
}
```

$P_0 = 0$ 意味着**「我对初始状态 100% 确定」**。这是个激进的选择 —— 严格说会让最初几步的卡尔曼增益 $K \approx 0$（完全不信测量）。好在每步预测都会 `+ Q`，$P$ 会迅速增长到合理量级，所以实际影响是**开头几步响应偏慢**。

对 `DKalmanObserver` 而言，$X_0 = [p(0);\, 0]$ —— 动量部分用实测值（确实准），扰动部分假设为 0（未必准，但 $Q$ 会很快把不确定性注入进来）。

#### ⚠️ 两个 `step()` 重载对 $A$ 的要求相反

这是**读这个类最容易踩的坑**：


| 重载             | 对`A` 的假设         | 内部处理                    |
| ---------------- | -------------------- | --------------------------- |
| `step(u, y)`     | `A` 已经是**离散**的 | 直接用                      |
| `step(u, y, dt)` | `A` 是**连续**的     | `At = I + dt*A`（前向欧拉） |

```cpp
Eigen::VectorXd KalmanFilter::step(Eigen::VectorXd& u, Eigen::VectorXd& y, double dt)
{
  At = I + dt*A;                           // ← 前向欧拉离散化
  X = At * X + dt * (B * u);
  P = At * P * At.transpose() + Q;         // ← ★ Q 没有乘 dt！
  ...
}
```

两个调用方都自洽，但方向相反：

- [tests/test_kalman.cpp:24](tests/test_kalman.cpp#L24) 手动离散化 `A = I + DELTA*A`，然后用 `step(u,y)` ✓
- [disturbance_kalman_filter.h:106](lib/disturbance_kalman_filter.h#L106) 传连续 `A`，用 `step(u,y,dt)` ✓

**接口本身没有任何提示** —— 两个重载名字一样，签名只差一个 `dt`。用错了不报错，只静默给出错误结果。

#### ⚠️ `Q` 未按 `dt` 缩放

上面代码里，`A` 和 `B` 都被 `dt` 正确缩放了（`I + dt*A`、`dt*(B*u)`），但：

```cpp
P = At * P * At.transpose() + Q;      // Q 原样使用
```

若 `Q` 代表连续时间过程噪声的**功率谱密度**，正确的离散化应该是 $Q_d \approx Q\,\Delta t$。这里没有除/乘 `dt`。

这不算 bug（`Q` 本来就是调参旋钮，可以理解成"每步的噪声协方差"），但有个实际后果：

> **`Q` 的含义是「每步」而不是「每秒」，所以控制周期一改，`Q` 必须重新整定。**

对照 `KalmanFilterContinous` 就能看出差异 —— 它用 Van Loan 正确积分了 $Q_d$，还做了 `Rd = Rupd / dt`。两个类对协方差的处理哲学不一致。

---

### 5.5 [lib/kalman_filter_continous.cpp](lib/kalman_filter_continous.cpp) —— Van Loan 精确离散化

`KalmanFilter` 用前向欧拉离散（$I + \Delta t A$），$\Delta t$ 稍大就不准。这个类改用**矩阵指数**做精确离散化。

#### 离散化 $A$ 和 $B$

连续系统 $\dot X = AX + Bu$ 的精确离散化是：

$$
A_d = e^{A\Delta t}, \qquad B_d = \int_0^{\Delta t} e^{A\tau}\,d\tau \; B
$$

$B_d$ 那个积分不好算。**技巧**：构造增广矩阵，一次矩阵指数同时得到两者：

$$
\exp\left(\begin{bmatrix} A & B \\ 0 & 0\end{bmatrix}\Delta t\right)
= \begin{bmatrix} e^{A\Delta t} & \int_0^{\Delta t}e^{A\tau}d\tau\,B \\ 0 & I \end{bmatrix}
= \begin{bmatrix} A_d & B_d \\ 0 & I \end{bmatrix}
$$

对应代码（构造函数里搭好 `AB`，`makeDiscrete` 里取块）：

```cpp
// 构造函数
AB.resize(na+nb, na+nb);
AB.setZero();                        // ← 清零
AB.block(0, 0, na, na) = a;          // 左上 = A
AB.block(0, na, na, nb) = b;         // 右上 = B
ABd.resize(na+nb, na+nb);

// makeDiscrete()
ABd = exponential(AB, dt);
Ad = ABd.block(0, 0, na, na);        // 左上 = A_d
Bd = ABd.block(0, na, na, nb);       // 右上 = B_d
```

#### 离散化 $Q$ —— Van Loan 方法

过程噪声协方差的精确离散化是个积分：

$$
Q_d = \int_0^{\Delta t} e^{A\tau}\,Q\,e^{A^T\tau}\,d\tau
$$

**Van Loan (1978)** 的技巧：构造另一个增广矩阵

$$
\mathcal{M} = \begin{bmatrix} -A & Q \\ 0 & A^T \end{bmatrix}\Delta t
$$

则

$$
e^{\mathcal{M}} = \begin{bmatrix} X_{11} & X_{12} \\ 0 & X_{22}\end{bmatrix}
= \begin{bmatrix} e^{-A\Delta t} & A_d^{-1}Q_d \\ 0 & A_d^T\end{bmatrix}
$$

于是

$$
\boxed{\;Q_d = X_{22}^T\,X_{12} = A_d\cdot\left(A_d^{-1}Q_d\right) = Q_d\;}\quad\checkmark
$$

**只用一次矩阵指数就拿到了积分结果，且不需要显式求 $A_d^{-1}$。**

对应代码（一行顶一篇论文）：

```cpp
AQd = exponential(AQ, dt);
Qd = AQd.block(na, na, na, na).transpose() * AQd.block(0, na, na, na);
//   └────── X₂₂ᵀ = A_d ──────┘            └───── X₁₂ = A_d⁻¹Q_d ─────┘
```

> ⚠️ **这里有个真实的 bug**：`AQ` 的左下块从未初始化，见 §10.1。上面推导中「$e^{\mathcal M}$ 是块上三角」这个性质**依赖左下块为零**，垃圾值会让整个结果失效。

#### `setCovariance` 的隐藏行为

```cpp
void KalmanFilterContinous::setCovariance(Eigen::MatrixXd& q, Eigen::MatrixXd& r)
{
  AQ.block(0, na, na, na) = q;     // ← Q 直接写进增广矩阵，没有单独存
  R = r;
}
```

注意 **`Q` 没有独立的成员变量** —— 它被直接嵌进 `AQ` 的右上块。而 `Qd`（离散化后的）才是成员。读代码时找不到 `Q` 成员不要疑惑。

#### 矩阵指数：5 项泰勒截断

```cpp
#define EXP_TERMS 5

Eigen::MatrixXd KalmanFilterContinous::exponential(Eigen::MatrixXd& m, double dt)
{
  Eigen::MatrixXd res = Eigen::MatrixXd::Identity(m.rows(), m.cols());
  Eigen::MatrixXd acc = res;
  for(int i = 1; i <= EXP_TERMS; i++) {
    acc *= m * (dt/i);       // acc = (m·dt)^i / i!  ← 递推，避免重复算幂和阶乘
    res += acc;
  }
  return res;
}
```

计算的是截断泰勒级数：

$$
e^{M\Delta t} \approx \sum_{i=0}^{5} \frac{(M\Delta t)^i}{i!}
$$

**递推很巧妙**：`acc *= m * (dt/i)` 让第 $i$ 项由第 $i{-}1$ 项乘上 $\frac{M\Delta t}{i}$ 得到，避免了重复计算矩阵幂和阶乘。

**但这不是精确的矩阵指数。** 截断误差约为

$$
\left\|\frac{(M\Delta t)^6}{6!}\right\| \approx \frac{\|M\Delta t\|^6}{720}
$$

注释掉的 `(AB*dt).exp()` 是 Eigen 的 `unsupported/MatrixFunctions` 模块（用 Padé 近似 + scaling-and-squaring，精确得多）—— 作者选择手写以避免额外依赖。

> **隐藏的适用性限制**：对小 $\|M\Delta t\|$（比如 $\Delta t = 0.01$、$A$ 的特征值 $O(1)$）够用；**若 $A$ 的特征值大或 $\Delta t$ 大，5 项会明显不准**，而且不会有任何警告。

#### `updateR()` —— 位形相关的测量噪声（很细致的设计）

```cpp
void KalmanFilterContinous::updateR(Eigen::MatrixXd& m) { Rupd = m * R * m.transpose(); }
```

**为什么需要它**：卡尔曼滤波的观测量是 $y = p = M(q)\dot q$，但真正带噪的是**速度** $\dot q$。这一行做的是**协方差通过线性变换的传播**。下面从头推。

**通用引理（协方差如何随线性变换改变）**：设随机向量 $x$ 的协方差为 $\Sigma_x$，$A$ 是**确定性**矩阵，$y = Ax$，则

$$
\mathrm{Cov}(y) = A\,\Sigma_x\,A^T
$$

*证明*（三步，从协方差定义出发）。记 $\bar x = E[x]$，则 $\bar y = E[Ax] = A\bar x$，偏差 $y - \bar y = A(x - \bar x)$。代入定义：

$$
\mathrm{Cov}(y)
= E\big[(y-\bar y)(y-\bar y)^T\big]
= E\big[A(x-\bar x)(x-\bar x)^T A^T\big]
= A\,\underbrace{E\big[(x-\bar x)(x-\bar x)^T\big]}_{\Sigma_x}\,A^T
$$

$A$ 是常数所以能提到期望外；$\big(A(x-\bar x)\big)^T = (x-\bar x)^T A^T$ 转置反序，这就是 $A^T$ 落在右边的原因。**此结论不需要高斯假设**，对任意分布的线性变换都精确成立。

**代入本问题**。速度测量带噪：

$$
\dot q_{meas} = \dot q + v, \qquad v \sim (0,\ R)
$$

$v$ 零均值、协方差 $R$（即用户给的、**速度层面**的测量噪声）。滤波器吃进去的"测量"是动量：

$$
p_{meas} = M(q)\,\dot q_{meas} = \underbrace{M(q)\dot q}_{p_{true}} + \underbrace{M(q)\,v}_{p_{noise}}
$$

噪声部分 $p_{noise} = M(q)\,v$ 正是「$v$ 经线性映射 $A = M(q)$」。套引理（$x\to v,\ \Sigma_x\to R,\ A\to M(q)$）：

$$
\boxed{\;\mathrm{Cov}(p_{noise}) = M(q)\,R\,M(q)^T\;}
$$

逐字对应代码 `Rupd = m * R * m.transpose()`：`m` = $M(q)$、`R` = 速度噪声、`Rupd` = 动量层面的噪声。

**测量噪声协方差随机器人位形变化。** [disturbance_kalman_filter_exp.h:110](lib/disturbance_kalman_filter_exp.h#L110) 每步调用 `filter->updateR(M)` 把这个效应喂进去。这是 `DKalmanObserverExp` 相对 `DKalmanObserver` 的**主要理论改进** —— 基础版偷懒用一个**常数** $R$ 假装它不变，但同样的速度噪声在惯量大的位形（手臂伸直）会被 $M$ 放大得多。

> **两个隐含假设**（读代码时要知道）：
> 1. **$M(q)$ 被当成确定性矩阵**。其实 $q$ 也有噪声、$M$ 也是随机的，但位置由编码器测、精度远高于（常需差分得到的）速度，所以近似 $M(q)$ 为已知常量，忽略 $q$ 噪声的二阶效应。
> 2. **$v$ 零均值**。否则 $p_{noise}$ 有偏，仅靠协方差不足以描述，还需单独处理偏置。

`makeDiscrete` 里还有一步：

```cpp
Rd = Rupd / dt;
```

连续时间测量噪声的**功率谱密度** $\to$ 离散采样协方差要除以采样周期：

$$
R_d = \frac{R_c}{\Delta t}
$$

标准做法。（对比 §5.4 —— `KalmanFilter` 完全没做这类缩放。）

---

## 6. 七个观测器逐个详解

每个小节的结构：**理论推导 → 成员变量表 → 逐行代码 → 首次调用路径 → 调参 → RNEA 版差异**。全部回到 §2.2 的核心方程 $\dot p = \tau - \beta + \tau_{ext}$。

---

### 6.1 [lib/momentum_observer.h](lib/momentum_observer.h) —— 动量观测器（De Luca）

**最经典、最常用、最便宜的一个。先读它。**

#### 理论

构造残差 $r$：

$$
r(t) = K_O\left[\;p(t) - \int_0^t \big(\tau - \beta + r\big)\,ds \;-\; p(0)\;\right]
$$

求导，代入核心方程 $\dot p = \tau - \beta + \tau_{ext}$：

$$
\dot r = K_O\left[\dot p - (\tau - \beta + r)\right]
       = K_O\left[(\tau - \beta + \tau_{ext}) - \tau + \beta - r\right]
$$

$$
\boxed{\;\dot r = K_O\,(\tau_{ext} - r)\;}
$$

**结论**：$r$ 是 $\tau_{ext}$ 的一阶低通滤波。逐关节的传递函数：

$$
\frac{r_i(s)}{\tau_{ext,i}(s)} = \frac{k_i}{s + k_i}
$$

时间常数 $1/k_i$。$k_i$ 越大响应越快，但对噪声和模型误差越敏感。

> **这个构造的精妙之处**：它完全不需要微分 $p$ —— 恰恰相反，把 $p$ 放在积分号**外面**，需要微分的量被积分抵消了。这是全库最优雅的想法。

#### 成员变量表

```cpp
private:
  VectorJ sum, r;                    // 累加器
  VectorJ p, beta, torque, tprev;    // 中间变量
  VectorJ ko;                        // 增益
```


| 变量     | 数学                       | 说明                                                |
| -------- | -------------------------- | --------------------------------------------------- |
| `p`      | $p = M\dot q$              | **会被复用**：后面变成 $p - \texttt{sum}$           |
| `beta`   | $\beta = G - C^T\dot q$    | 核心方程里的组合项                                  |
| `sum`    | $p(0) + \int_0^t(\cdot)ds$ | **身兼二职**：积分累加器 **且** 初始动量            |
| `r`      | $r$                        | 残差，即$\tau_{ext}$ 的估计                         |
| `torque` | 多义                       | **复用三次**：$\tau - \tau_f$ → 被积函数 → 返回值 |
| `tprev`  | 上一步的被积函数           | 梯形积分需要                                        |
| `ko`     | $K_O$                      | 逐关节增益（向量，非矩阵）                          |

#### 逐行代码

```cpp
VectorJ MomentumObserver::getExternalTorque(VectorJ& q, VectorJ& qd, VectorJ& tau, double dt)
{
  p = dyn->getM(q) * qd;                                    // ① p = M·q̇
  beta = dyn->getG(q) - dyn->getC(q, qd).transpose() * qd;  // ② β = G - Cᵀq̇

  torque = tau - dyn->getFriction(qd);                      // ③ τ - τ_f

  if(isRun) {
    torque += r - beta;                    // ④ 被积函数 = τ - τ_f + r - β
    sum += 0.5 * dt * (torque + tprev);    // ⑤ 梯形积分累加
  } else {
    torque -= beta;                        // ⑥ 首次：被积函数（r=0）
    r.setZero();
    sum = p;                               // ⑦ ★ sum 初始化为 p(0)
    isRun = true;
  }
  tprev = torque;                          // ⑧ 保存本步被积函数

  p -= sum;                                // ⑨ p - p(0) - ∫

  for(int i = 0; i < jointNo; i++) {
    r(i) = ko(i) * p(i);                   // ⑩ r = K_O ⊙ (...)  逐元素积
    torque(i) = r(i);                      // ⑪ 复用 torque 当返回值
  }
  return torque;
}
```

**关键行解释**：

**④ 被积函数为什么是 $\tau - \tau_f + r - \beta$？**
理论公式里积分号内是 $\tau - \beta + r$。代码里的 `tau` 还含摩擦，所以先扣掉：$\tau_{\text{真实驱动}} = \tau_{\text{测量}} - \tau_f$。

**⑤ 梯形积分**：

$$
\int_{t_k}^{t_{k+1}} f\,ds \approx \frac{\Delta t}{2}\big(f(t_k) + f(t_{k+1})\big)
$$

精度 $O(\Delta t^2)$，比前向欧拉的 $O(\Delta t)$ 高一阶。**对积分器类观测器很重要** —— 误差会永久累积，不像滤波器那样自动衰减。

**⑦⑨ `sum` 身兼二职（★最容易看晕的地方）**：
首次调用 `sum = p`，即 $\texttt{sum} = p(0)$。之后每步 `sum += 梯形增量`。所以任意时刻：

$$
\texttt{sum} = p(0) + \int_0^t(\tau - \beta + r)\,ds
$$

于是第 ⑨ 行 `p -= sum` 一次性完成了 $p(t) - p(0) - \int$，正好是理论公式方括号里的内容。

省了一个变量，但让 `sum` 这个名字很误导 —— 它不只是「和」。

**⑩ 逐元素积不是矩阵乘**：`ko(i) * p(i)` 说明 $K_O$ 是**对角**增益矩阵，用向量存储。每个关节独立整定。

#### 首次调用路径追踪

```
① p    = M·q̇  = p(0)
② beta = β(0)
③ torque = τ(0) - τ_f(0)
⑥ torque -= beta      → τ(0) - τ_f(0) - β(0)
   r.setZero()        → r = 0
⑦ sum = p             → sum = p(0)
⑧ tprev = torque      → 保存，供第二步的梯形积分用
⑨ p -= sum            → p = p(0) - p(0) = 0     ★
⑩ r(i) = ko(i)·0 = 0
⑪ torque(i) = 0
   return 0  ✓
```

**第二次调用**验证梯形积分的正确性：

```
④ torque = τ(1) - τ_f(1) + r(0) - β(1)  = f(1)     (r(0)=0)
⑤ sum = p(0) + 0.5·dt·(f(1) + f(0))                ✓ 标准梯形
⑨ p = p(1) - p(0) - 0.5·dt·(f(1)+f(0))             ✓
```

#### 调参

```cpp
VectorJ k(6);  k << 30,30,30,30,30,30;
MomentumObserver obs(&robot, k);
```

`k` 是逐关节增益，单位是 $\mathrm{s}^{-1}$（时间常数的倒数）。经验值 $30\sim50$：


| $k$ | 时间常数$1/k$ | 效果                              |
| --- | ------------- | --------------------------------- |
| 10  | 100 ms        | 太慢，碰撞检测延迟大              |
| 30  | 33 ms         | `test.cpp` 用这个                 |
| 50  | 20 ms         | `test.m` / `test_rnea.cpp` 用这个 |
| 200 | 5 ms          | 太快，放大模型误差和噪声 → 误报  |

**整定思路**：$k$ 是「响应速度」和「抗噪声/抗模型误差」的直接权衡。从 30 开始，观察静止时的残差噪声水平，在噪声可接受的前提下往大调。

#### RNEA 版差异

[momentum_observer_rnea.h:75-76](lib/momentum_observer_rnea.h#L75-L76) 只有两行不同：

```cpp
// 原版
p    = dyn->getM(q) * qd;
beta = dyn->getG(q) - dyn->getC(q, qd).transpose() * qd;

// RNEA 版
p    = dyn->rnea(q, zero, qd);                                   // M·q̇
beta = dyn->rnea(q, zero, zero, GRAVITY) - dyn->tranCqd(q, qd);  // G - Cᵀq̇
```

其余逻辑**逐字相同**。多了一个 `zero` 成员（`VectorJ::Zero(jointNo)`）作为常量零向量。

> **七个观测器的 RNEA 版本都是这个模式** —— 只替换动力学调用，控制逻辑不变。读懂这一对，其余四对扫一眼即可。

---

### 6.2 [lib/disturbance_observer.h](lib/disturbance_observer.h) —— 扰动观测器（Mohammadi）

> ⚠️ **勘误说明**：本节的 $\dot z$ 方程和收敛性证明在早先版本里符号是错的（$C\dot q$、$G$、$\tau$ 三项符号全反，导致"证明"最后一步无法成立）。下面是**重新对照代码逐步推导、并用稳态值交叉验证过**的正确版本。

#### 问题设定

先明确要估计什么。带外力的动力学（暂略摩擦，最后再加回）：

$$
M(q)\ddot q + C(q,\dot q)\dot q + G(q) = \tau + \tau_{ext}
$$

于是外部力矩

$$
\tau_{ext} = M\ddot q + C\dot q + G - \tau \tag{$*$}
$$

反解出加速度项，这一步是后面代入的关键：

$$
M\ddot q = \tau + \tau_{ext} - C\dot q - G \tag{$**$}
$$

> **注意：这里出现的是 $C\dot q$，不是 $C^T\dot q$。** 和动量观测器不同 —— 那边因为要微分 $p = M\dot q$ 产生了 $\dot M$，才用 $\dot M = C + C^T$ 的恒等式引出转置。**这个观测器的辅助量是 $p = Y\dot q$（$Y$ 常数），$\dot p = Y\ddot q$ 里根本不出现 $\dot M$，所以全程用普通 $C$**。代码里也确实是 `getC(q,qd)*qd`（不转置），与动量观测器的 `getC(q,qd).transpose()*qd` 形成对照。

#### 辅助变量与设计自由度

Mohammadi 扰动观测器的核心思想：引入辅助向量 $p(\dot q)$ 和内部状态 $z$，令估计值

$$
\hat\tau_{ext} = z + p(\dot q)
$$

**为什么要拆成 $z + p$？** 因为直接估 $\tau_{ext}$ 需要 $\ddot q$（见 $(*)$），而 $\ddot q$ 测不到。把「含 $\ddot q$ 的部分」用 $p(\dot q)$ 的时间导数隐式吸收掉，$z$ 只需一阶积分 —— 又是那个「规避微分」的老套路。

关键的**设计关系**：$p$ 只是 $\dot q$ 的函数，所以链式法则

$$
\dot p = \frac{\partial p}{\partial \dot q}\,\ddot q
$$

代入 $(**)$，让待定的增益矩阵 $L$ 恰好吸收掉 $M^{-1}$：

$$
L(q) \triangleq \frac{\partial p}{\partial \dot q}\,M^{-1}
$$

代码选了**最简单的线性辅助量** $p = Y\dot q$（$Y$ 为常数矩阵），于是 $\dfrac{\partial p}{\partial \dot q} = Y$：

$$
\boxed{\;L = Y M^{-1}\;}
$$

正是代码第一行 `L = Y * dyn->getM(q).inverse()`。$Y$ 是自由选的常数增益（代码里 $Y = kI$）。

#### 推导 $\dot z$（逐步，符号已核对）

目标：设计 $\dot z$，使误差 $e = \tau_{ext} - \hat\tau_{ext}$ 满足漂亮的一阶动态 $\dot e = -Le + \dot\tau_{ext}$。

先算 $\dot p$。由 $p = Y\dot q$ 和 $(**)$：

$$
\dot p = Y\ddot q = \underbrace{YM^{-1}}_{L}\,M\ddot q = L\big(\tau + \tau_{ext} - C\dot q - G\big)
$$

现在要求

$$
\frac{d\hat\tau_{ext}}{dt} = L\big(\tau_{ext} - \hat\tau_{ext}\big)
$$

左边展开 $\dfrac{d\hat\tau_{ext}}{dt} = \dot z + \dot p$，代入上面的 $\dot p$，并令其等于目标：

$$
\dot z + L\big(\tau + \tau_{ext} - C\dot q - G\big) = L\tau_{ext} - L\underbrace{(z + p)}_{\hat\tau_{ext}}
$$

解出 $\dot z$（$L\tau_{ext}$ 两边抵消）：

$$
\dot z = -Lz - Lp - L\tau + LC\dot q + LG
$$

$$
\boxed{\;\dot z = -Lz - L\big(p + \tau - C\dot q - G\big)\;}
\qquad\text{等价于}\qquad
\dot z = -Lz + L\big(C\dot q + G - \tau - p\big)
$$

**右边这个形式和代码逐字对应** —— 代码的 `rht = z + L*(C·qd + G - τ - p)` 就是它的后向欧拉离散（见下）。

#### 收敛性证明（完整）

把上面设计的 $\dot z$ 代回 $\dfrac{d\hat\tau_{ext}}{dt} = \dot z + \dot p$：

$$
\begin{aligned}
\frac{d\hat\tau_{ext}}{dt}
&= \underbrace{\big[-Lz - L(p + \tau - C\dot q - G)\big]}_{\dot z}
 + \underbrace{L\big(\tau + \tau_{ext} - C\dot q - G\big)}_{\dot p} \\[4pt]
&= -Lz - Lp - L\tau + LC\dot q + LG + L\tau + L\tau_{ext} - LC\dot q - LG \\[4pt]
&= -L(z + p) + L\tau_{ext} \\[4pt]
&= -L\hat\tau_{ext} + L\tau_{ext}
\end{aligned}
$$

$$
\boxed{\;\frac{d\hat\tau_{ext}}{dt} = L\big(\tau_{ext} - \hat\tau_{ext}\big)\;}
$$

（这一次 $C\dot q$、$G$、$\tau$ 项**真的**两两抵消了 —— 对比早先错误版本硬凑的最后一步。）

误差 $e = \tau_{ext} - \hat\tau_{ext}$ 满足：

$$
\dot e = \dot\tau_{ext} - L\,e \quad\Longrightarrow\quad \dot e + Le = \dot\tau_{ext}
$$

若外力缓变（$\dot\tau_{ext} \approx 0$），则 $\dot e = -Le$，指数收敛。收敛速率由 $L = kM^{-1}$ 的特征值决定：$M$ 对称正定、$k>0$，故 $kM^{-1}$ 也对称正定，特征值全为正 —— **稳定 ✓**。

形态和动量观测器一样（一阶收敛），**区别是收敛速率矩阵 $L = YM^{-1}$ 随位形变化** —— 惯量大的位形收敛慢，惯量小的位形收敛快。这既是特点也是麻烦（性能不一致）。

#### 稳态交叉验证（确认符号无误）

令 $\dot z = 0$，从 $\dot z = -Lz + L(C\dot q + G - \tau - p)$ 得（$L$ 可逆）：

$$
z^{*} = C\dot q + G - \tau - p
$$

于是估计值

$$
\hat\tau_{ext} = z^{*} + p = C\dot q + G - \tau
$$

而由 $(*)$ 在稳态（$\ddot q = 0$）下 $\tau_{ext} = C\dot q + G - \tau$。**两者一致 ✓** —— 这就是符号正确的独立佐证。（用早先的错误公式做同样的稳态分析会得到 $\hat\tau_{ext} = -(C\dot q + G) - \tau$ 之类对不上的结果。）

#### 隐式离散化

对 $\dot z = -Lz + L(C\dot q + G - \tau - p)$ 用**后向欧拉**（$L$ 这里指 $L_0 = YM^{-1}$，尚未乘 $\Delta t$）：

$$
z_{k+1} = z_k + \Delta t\Big[-L_0\,z_{k+1} + L_0\big(C\dot q + G - \tau - p\big)\Big]
$$

把含 $z_{k+1}$ 的项归并到左边：

$$
\boxed{\;(I + \Delta t\,L_0)\,z_{k+1} = z_k + \Delta t\,L_0\big(C\dot q + G - \tau - p\big)\;}
$$

记 $\bar L = \Delta t\,L_0$（即代码里 `L *= dt` 之后的 `L`），就是

$$
(I + \bar L)\,z_{k+1} = z_k + \bar L\big(C\dot q + G - \tau - p\big)
$$

三块一一对应代码：`lft = I + L`、`rht = z + L*(C·qd + G - torque - p)`、`z = lft.inverse()*rht`。

**必须解一个线性方程组** —— 这就是代码里那次矩阵求逆的来源。

#### 成员变量表

```cpp
private:
  MatrixJ Y, L, I;
  MatrixJ lft, rht;
  VectorJ p, z, torque;
```


| 变量  | 数学                          | 说明                                            |
| ----- | ----------------------------- | ----------------------------------------------- |
| `Y`   | $Y = k\,I$                    | 常数增益矩阵，由`settings()` 算出               |
| `L`   | $\Delta t\cdot YM^{-1}$       | **已吸收 $\Delta t$**（`L *= dt`）              |
| `I`   | $I$                           | 单位阵，构造时算好                              |
| `lft` | $I + \Delta t L$              | 线性方程组左边                                  |
| `rht` | 右边                          | 注意它声明为`MatrixJ` 但存的是向量（见下）      |
| `p`   | $Y\dot q$ → $\hat\tau_{ext}$ | **复用**：先是 $Y\dot q$，最后加上 $z$ 变返回值 |
| `z`   | $z$                           | 观测器内部状态                                  |

#### 逐行代码

```cpp
VectorJ DisturbanceObserver::getExternalTorque(VectorJ& q, VectorJ& qd, VectorJ& tau, double dt)
{
  L = Y * dyn->getM(q).inverse();          // ① L = Y·M⁻¹        ← 第 1 次求逆
  L *= dt;                                 // ② L ← Δt·L          ★ 吸收 dt
  p = Y * qd;                              // ③ p = Y·q̇

  if(isRun) {
    torque = tau - dyn->getFriction(qd);   // ④ τ - τ_f
    lft = I + L;                           // ⑤ I + Δt·L
    rht = z + L*(dyn->getC(q, qd)*qd + dyn->getG(q) - torque - p);
    //     └z_k┘  └Δt·L┘ └── C·q̇ + G - (τ-τ_f) - Y·q̇ ──┘        ⑥
    z = lft.inverse() * rht;               // ⑦ 解方程组          ← 第 2 次求逆
  } else {
    z = -p;                                // ⑧ 使 τ̂_ext(0) = 0
    isRun = true;
  }

  p += z;                                  // ⑨ τ̂_ext = Y·q̇ + z
  return p;
}
```

**关键行解释**：

**② `L *= dt` —— 全函数最容易看漏的一行。**
之后所有出现 `L` 的地方都已经含 $\Delta t$。所以第 ⑤ 行 `lft = I + L` 就是 $I + \Delta t L$，第 ⑥ 行 `L*(...)` 就是 $\Delta t L(\cdots)$。**如果没注意到这一行，整个离散化会看不懂。**

**⑥ 对照理论公式**：
上一节推出的离散方程右边是 $z_k + \bar L\big(C\dot q + G - \tau - p\big)$。代码是 `z + L*(C*qd + G - torque - p)`，其中 `torque = τ - τ_f`、`p = Y·qd`、`L = `$\bar L$。代入 `torque`：

$$
z + \bar L\big(C\dot q + G - (\tau - \tau_f) - Y\dot q\big)
= z + \bar L\big(\underbrace{C\dot q + G - \tau - p}_{\text{无摩擦时的项}} + \tau_f\big)
$$

**逐字对上 ✓**。带摩擦时只是把 $\tau$ 换成 $\tau - \tau_f$（先扣掉摩擦再当作驱动力矩），符号方向和 $(*)$、$(**)$ 完全一致。

**⑦ 两次矩阵求逆**：第 ① 行的 $M^{-1}$ 和第 ⑦ 行的 $(I+\Delta t L)^{-1}$。**这是全库除 `DKalmanObserverExp` 外最贵的观测器。**

> 数值上更好的写法是 `z = lft.lu().solve(rht);`（LU 分解解方程组）而不是显式求逆 —— 更快也更稳。对 $6\times6$ 差别不大，但这是个可改进点。

**⑧ 首次初始化**：$z = -Y\dot q$，于是第 ⑨ 行 $\hat\tau_{ext} = Y\dot q + (-Y\dot q) = 0$ ✓。

**隐式 vs 显式**：


|                      |                                                   |
| -------------------- | ------------------------------------------------- |
| ✅**无条件稳定**     | $\Delta t$ 再大也不发散 —— 这是选隐式的全部理由 |
| ❌ 每步 2 次矩阵求逆 | $O(n^3)$ × 2                                     |

对比 `SlidingModeObserver` 用的是**显式**前向欧拉，$\Delta t$ 太大会崩。

#### 首次调用路径追踪

```
① L = Y·M⁻¹        ← 算了但没用上（浪费一次矩阵求逆）
② L *= dt          ← 同样浪费
③ p = Y·q̇
⑧ z = -p           → z = -Y·q̇
⑨ p += z           → p = Y·q̇ - Y·q̇ = 0      ★
   return 0  ✓
```

> **首次调用白算了一次 $M^{-1}$** —— `L` 在 `else` 分支里完全没用到。无害，但是个小浪费。

#### 调参 —— 三个参数塌缩成一个

```cpp
void DisturbanceObserver::settings(double sigma, double xeta, double beta)
{
  double k = 0.5*(xeta + 2*beta*sigma);
  Y = k * MatrixJ::Identity(jointNo, jointNo);
}
```

$$
k = \frac{1}{2}\big(\xi + 2\beta\sigma\big), \qquad Y = k\,I
$$

**三个参数最后只塌缩成一个标量** $k$。三个自由度实际只有一个 —— 无数组合给出同一个 $Y$。

它们的 Doxygen 注释全是 `@param sigma ...` —— **作者没写**。要理解只能翻 Mohammadi 的论文。`test.cpp` / `test.m` 里的经验值：

```cpp
double sigma = 21, xeta = 18, beta = 50;
// → k = 0.5*(18 + 2*50*21) = 0.5*2118 = 1059
```

> **实用建议**：既然只有 $k$ 起作用，调参时固定 `sigma=1, xeta=0`，只调 `beta`（此时 $k = \beta$），比调三个耦合参数直观得多。

#### RNEA 版差异

[disturbance_observer_rnea.h:86](lib/disturbance_observer_rnea.h#L86)：

```cpp
// 原版
rht = z + L*(dyn->getC(q, qd)*qd + dyn->getG(q) - torque - p);
// RNEA 版
rht = z + L*(dyn->rnea(q, qd, zeros, GRAVITY) - torque - p);
//            └── 一次调用同时拿到 C·q̇ + G ──┘
```

`getM()` 仍然要调用（基类实现，$N$ 次 RNEA）。这是唯一一个用 `getM()` 而非 `tranCqd()` 的 RNEA 观测器，所以只要 7 次调用而不是 10 次。

---

### 6.3 [lib/sliding_mode_observer.h](lib/sliding_mode_observer.h) —— 滑模观测器（Garofalo 2019）

理论上最强（**有限时间**收敛，不是渐近），但参数最多、最容易抖振。

> 论文：Garofalo et al., *Sliding mode momentum observers for estimation of external torques and joint acceleration*, 2019.

#### 理论 —— super-twisting 二阶滑模

估计动量 $\hat p$，定义误差：

$$
\tilde p = \hat p - p
$$

观测器动态：

$$
\begin{aligned}
\dot{\hat p} &= \tau + C^T\dot q - G + \sigma - T_2\,\tilde p - T_1\sqrt{|\tilde p|}\,\mathrm{sign}(\tilde p) \\
\dot\sigma &= -S_1\,\mathrm{sign}(\tilde p) - S_2\,\tilde p
\end{aligned}
$$

**误差动态**（代入真实的 $\dot p = \tau - \beta + \tau_{ext} = \tau + C^T\dot q - G + \tau_{ext}$）：

$$
\begin{aligned}
\dot{\tilde p} &= \dot{\hat p} - \dot p = \sigma - T_2\tilde p - T_1\sqrt{|\tilde p|}\,\mathrm{sign}(\tilde p) - \tau_{ext} \\
\dot\sigma &= -S_1\mathrm{sign}(\tilde p) - S_2\tilde p
\end{aligned}
$$

这是 **super-twisting 算法**的标准结构。它保证 $\tilde p \to 0$ 且 $\sigma \to \tau_{ext}$，**在有限时间内**（不是渐近）。

$$
\boxed{\;\sigma \text{ 就是外部力矩的估计，代码返回它}\;}
$$

**为什么 $\sigma \to \tau_{ext}$**：滑模到达后 $\tilde p \equiv 0$，$\dot{\tilde p} \equiv 0$，代入第一式得 $\sigma - \tau_{ext} = 0$。

**线性项的作用**：$-T_2\tilde p$ 和 $-S_2\tilde p$ 是**额外的**线性项（"generalized super-twisting"），可以设 $S_2 = T_2 = 0$ 退化成经典 super-twisting。[test.cpp:42](tests/test.cpp#L42) 的注释 `// set 0 to exclude linear part` 就是说这个。加上线性项能改善远离滑模面时的收敛速度。

#### 成员变量表

```cpp
private:
  const double BIG = 50.0;             /**< Map tanh to sign. */
  VectorJ sigma, p_hat;                // 观测器状态（持久）
  VectorJ p, spp, dp_hat, torque, dsigma;   // 每步临时量
  VectorJ T1, S1, T2, S2;              // 增益
```


| 变量          | 数学                      | 说明                                      |
| ------------- | ------------------------- | ----------------------------------------- |
| `sigma`       | $\sigma$                  | **状态**，也是最终答案                    |
| `p_hat`       | $\hat p$                  | **状态**，动量估计                        |
| `p`           | 多义                      | **复用三次**：$p$ → $\tilde p$ → 返回值 |
| `spp`         | $\mathrm{sign}(\tilde p)$ | 用$\tanh$ 平滑逼近                        |
| `dp_hat`      | $\dot{\hat p}$            |                                           |
| `dsigma`      | $\dot\sigma$              |                                           |
| `T1,S1,T2,S2` |                           | 逐关节增益向量                            |

#### 逐行代码

```cpp
VectorJ SlidingModeObserver::getExternalTorque(VectorJ& q, VectorJ& qd, VectorJ& tau, double dt)
{
  p = dyn->getM(q) * qd;                 // ① p = M·q̇
  torque = tau - dyn->getFriction(qd);   // ② τ - τ_f

  if(!isRun) {
    sigma.setZero();                     // ③ σ(0) = 0
    p_hat = p;                           // ④ p̂(0) = p(0)  → p̃(0) = 0
    isRun = true;
  }

  p = p_hat - p;                         // ⑤ ★ p 复用为误差 p̃
  for(int i = 0; i < jointNo; i++) {
    spp(i) = tanh(p(i)*BIG);             // ⑥ ★ tanh(50·p̃) ≈ sign(p̃)
  }

  dp_hat = torque + dyn->getC(q, qd).transpose()*qd - dyn->getG(q) + sigma;
  //       └── τ - τ_f + Cᵀq̇ - G + σ ──┘                                  ⑦
  for(int i = 0; i < jointNo; i++) {
    dp_hat(i) -= T2(i)*p(i);                          // ⑧ -T₂·p̃
    dp_hat(i) -= sqrt(fabs(p(i))) * T1(i) * spp(i);   // ⑨ -T₁·√|p̃|·sign(p̃)
  }
  for(int i = 0; i < jointNo; i++) {
    dsigma(i) = -S1(i)*spp(i) - S2(i)*p(i);           // ⑩ σ̇
  }

  p = sigma;                             // ⑪ ★ p 再次复用，保存返回值
  p_hat += dt*dp_hat;                    // ⑫ 前向欧拉
  sigma += dt*dsigma;                    // ⑬ 前向欧拉
  return p;                              // ⑭ 返回更新前的 σ
}
```

**关键行解释**：

**⑤⑪ `p` 被复用三次（★全库可读性最差的地方）**：


| 行    | `p` 的含义                          |
| ----- | ----------------------------------- |
| ① 后 | 广义动量$p = M\dot q$               |
| ⑤ 后 | 动量估计误差$\tilde p = \hat p - p$ |
| ⑪ 后 | 返回值（外部力矩估计$\sigma$）      |

**同一个变量三种完全不同的物理量。** 读的时候必须盯紧 `// reuse` 注释。

**⑥ `tanh(50·x)` 平滑逼近 `sign(x)`**：

$$
\mathrm{sign}(x) \approx \tanh(\kappa x), \qquad \kappa = 50
$$

真正的 `sign()` 在 $\tilde p = 0$ 附近会引起**抖振**（chattering）—— 数值上高频切换，物理上激励机械共振、磨损减速器。用 $\tanh$ 抹平是标准工程做法。

被注释掉的那行就是真 `sign`：

```cpp
// for(int i = 0; i < jointNo; i++) spp(i) = (p(i) > 0 ? 1 : (p(i) < 0 ? -1 : 0));
```

$\kappa$ 的权衡：越大越接近真 sign（收敛更快、理论保证更强），但抖振越严重。$\kappa = 50$ 意味着过渡带宽约 $2/\kappa = 0.04$（$\tanh$ 在 $|x| < 0.04$ 内近似线性）。

**⑨ `sqrt(fabs(p(i)))`**：super-twisting 的标志性非线性项 $\sqrt{|\tilde p|}$。它是有限时间收敛的来源 —— 在 $\tilde p \to 0$ 时，$\sqrt{|\tilde p|}$ 的导数趋于无穷，产生"无限强"的收敛驱动力。

**⑫⑬ 前向欧拉（显式）**：

$$
\hat p_{k+1} = \hat p_k + \Delta t\,\dot{\hat p}_k, \qquad
\sigma_{k+1} = \sigma_k + \Delta t\,\dot\sigma_k
$$

**和 `DisturbanceObserver` 的隐式相反** —— $\Delta t$ 太大会不稳定。滑模系统对步长尤其敏感（因为高增益切换），实践中需要较高的控制频率。

**⑭ 返回更新前的 $\sigma$**：`p = sigma;`（第 ⑪ 行）在 `sigma += dt*dsigma;`（第 ⑬ 行）**之前**。所以返回值滞后一个时间步。无关紧要，但读代码时会疑惑「为什么不返回最新的」。

#### 首次调用路径追踪

```
① p = M·q̇ = p(0)
② torque = τ(0) - τ_f(0)
③ sigma = 0
④ p_hat = p = p(0)
⑤ p = p_hat - p = 0                    → p̃ = 0  ★
⑥ spp = tanh(0) = 0
⑦ dp_hat = torque + Cᵀq̇ - G + 0
⑧ dp_hat -= T2·0 = 0                   （无影响）
⑨ dp_hat -= sqrt(0)·T1·0 = 0           （无影响）
⑩ dsigma = -S1·0 - S2·0 = 0
⑪ p = sigma = 0                        ★
⑫ p_hat += dt·dp_hat
⑬ sigma += dt·0 = 0
   return 0  ✓
```

#### 调参 —— 4 个向量，实际调 2 个

```cpp
VectorJ T1(2), S1(2), T2(2), S2(2);
S1 << 15,20;
T1(0) = 2*sqrt(S1(0));  T1(1) = 2*sqrt(S1(1));    // ★ T₁ = 2√S₁
S2 << 10,10;                                       // 设 0 关掉线性部分
T2(0) = 2*sqrt(S2(0));  T2(1) = 2*sqrt(S2(1));    // ★ T₂ = 2√S₂
```

$$
T_1 = 2\sqrt{S_1}, \qquad T_2 = 2\sqrt{S_2}
$$

**这是 super-twisting 的经典整定规则**（保证收敛的充分条件之一）。所以实际只需要调：

- **$S_1$**：非线性部分强度 —— 决定有限时间收敛的速度，也决定对扰动的鲁棒性上界
- **$S_2$**：线性部分强度 —— 设 0 退化成经典 super-twisting

各处的经验值：


| 来源                       | $S_1$      | $S_2$      |
| -------------------------- | ---------- | ---------- |
| `test.cpp`                 | `[15, 20]` | `[10, 10]` |
| `test.m` / `test_rnea.cpp` | `[20, 30]` | `[10, 10]` |

#### RNEA 版差异

[sliding_mode_observer_rnea.h:94,108](lib/sliding_mode_observer_rnea.h#L94)：

```cpp
p = dyn->rnea(q, zeros, qd);
...
dp_hat = torque + dyn->tranCqd(q, qd) - dyn->rnea(q, zeros, zeros, GRAVITY) + sigma;
//                └── Cᵀq̇ ──┘         └────── G ──────┘
```

---

### 6.4 [lib/disturbance_kalman_filter.h](lib/disturbance_kalman_filter.h) —— 卡尔曼扰动观测器（Hu）

#### 理论 —— 把 $\tau_{ext}$ 建成状态

前三个观测器（动量、扰动、滑模）都是**确定性**的：设计一个误差动态让它收敛。这个观测器换了个哲学 —— 把外力估计变成**随机状态估计**问题，用卡尔曼滤波在过程噪声和测量噪声之间做统计最优的折中。

##### 第一步：还是那个核心方程

从 §2.2 的核心方程出发（这次把摩擦写进去，因为代码扣了摩擦）。带外力和摩擦的动力学，有效驱动力矩是 $\tau - \tau_f$：

$$
M\ddot q + C\dot q + G = (\tau - \tau_f) + \tau_{ext}
$$

对动量 $p = M\dot q$ 求导，代入上式并用斜对称恒等式 $\dot M = C + C^T$：

$$
\begin{aligned}
\dot p &= M\ddot q + \dot M\dot q
= \big[(\tau - \tau_f) + \tau_{ext} - C\dot q - G\big] + (C + C^T)\dot q \\
&= \underbrace{(\tau - \tau_f) - G + C^T\dot q}_{\displaystyle u\ (\text{全部已知})} + \tau_{ext}
\end{aligned}
$$

定义**已知输入** $u$，得到一个极简的一阶系统：

$$
\boxed{\;\dot p = u + \tau_{ext}\;}, \qquad u = (\tau - \tau_f) - G + C^T\dot q = \tau - \tau_f - \beta
$$

> 注意 $u = \tau - \tau_f - \beta$，$\beta = G - C^T\dot q$ 和动量观测器里完全同一个 $\beta$。代码 `u = tau - getG - getFriction; u += getC.transpose()*qd` 逐项对应 $-G$ 和 $+C^T\dot q$，符号正确。

这一步告诉我们：**动量的变化率 = 已知的 $u$ + 未知的外力**。如果能测到 $\dot p$，外力就是 $\dot p - u$ —— 但又回到微分噪声问题。卡尔曼的办法是不去微分，而是把 $\tau_{ext}$ 建成状态变量来估计。

##### 第二步：给外力配一个"外部模型"

关键的新思想（源自内模原理 / exosystem）：假设外力本身是某个**自治线性系统的输出**：

$$
\dot\omega = S\,\omega, \qquad \tau_{ext} = H\,\omega
$$

- $\omega \in \mathbb{R}^n$ 是扰动的**内部状态**
- $S$（`s` 参数）描述扰动**怎么随时间演化**
- $H$（`h` 参数）把内部状态映射成实际力矩

**为什么要多此一举？** 因为卡尔曼滤波需要一个「状态怎么演化」的模型才能做预测。直接说「$\tau_{ext}$ 是任意的」没法预测；说「$\tau_{ext}$ 服从 $\dot\omega = S\omega$」就给了滤波器一个先验。不同的 $S$ 编码不同的先验假设：

| $S$ | $H$ | 假设的外力形态 |
|---|---|---|
| $0$ | $I$ | **常值 / 随机游走**（最常用，测试代码全用它）：$\dot\omega = 0$，外力只由过程噪声驱动突变 |
| $\begin{bmatrix}0&1\\-\omega_0^2&0\end{bmatrix}$ | $[I\ 0]$ | 已知频率 $\omega_0$ 的**振动** |
| $-\lambda I$ | $I$ | 指数衰减的冲击 |

##### 第三步：拼成增广状态空间

把「物理」的 $\dot p = u + \tau_{ext} = u + H\omega$ 和「外部模型」的 $\dot\omega = S\omega$ 叠起来，状态取 $X = \begin{bmatrix}p\\\omega\end{bmatrix} \in \mathbb{R}^{2n}$：

$$
\underbrace{\begin{bmatrix}\dot p\\\dot\omega\end{bmatrix}}_{\dot X}
=
\underbrace{\begin{bmatrix}0 & H\\ 0 & S\end{bmatrix}}_{A}
\underbrace{\begin{bmatrix}p\\\omega\end{bmatrix}}_{X}
+
\underbrace{\begin{bmatrix}I\\ 0\end{bmatrix}}_{B} u,
\qquad
y = \underbrace{\begin{bmatrix}I & 0\end{bmatrix}}_{C} X = p
$$

逐块验证：
- 第一行 $\dot p = 0\cdot p + H\omega + Iu = u + H\omega$ ✓
- 第二行 $\dot\omega = 0\cdot p + S\omega + 0\cdot u = S\omega$ ✓
- 输出 $y = p$，即把动量当"测量" ✓

$$
A = \begin{bmatrix} 0 & H \\ 0 & S \end{bmatrix},\qquad
B = \begin{bmatrix} I \\ 0 \end{bmatrix},\qquad
C = \begin{bmatrix} I & 0 \end{bmatrix},\qquad
y = p = M(q)\dot q
$$

> **注意 $A$ 依赖 $q$ 吗？** 不 —— $A, B, C$ **全是常数矩阵**（$H, S$ 是用户给的常数）。位形相关性全被塞进了输入 $u$ 和测量 $y$ 里（它们每步用当前 $M, C, G$ 重算）。这正是为什么构造函数里能一次性建好 `KalmanFilter(A,B,C)`，之后只喂 $u$ 和 $y$。这是个漂亮的解耦。

##### 关于"测量" $y = p$ 的本质（重要）

$y = p = M(q)\dot q$ **不是直接测到的** —— 真正测的是 $q, \dot q$，$p$ 是乘 $M(q)$ 算出来的。这带来一个被 `DKalmanObserver` 忽略、被 `DKalmanObserverExp` 修正的问题：

真正带噪的是 $\dot q$（设协方差 $R_{\dot q}$）。经过 $M(q)$ 线性变换后，$p$ 的噪声协方差是**位形相关**的：

$$
\mathrm{Cov}(p_{noise}) = M(q)\,R_{\dot q}\,M(q)^T
$$

（这个协方差传播公式的完整推导 —— 通用引理及其证明、代入本问题、两个隐含假设 —— 见 §5.5 的 `updateR()` 小节。）

`DKalmanObserver` 用一个**常数** $R$ 近似它（假装测量噪声不随位形变化），`DKalmanObserverExp` 才用 §5.5 的 `updateR(M)` 每步算 $MRM^T$。这是两者的核心区别，见 §6.5。

#### 成员变量表

```cpp
private:
  MatrixJ H;                  // 扰动观测矩阵
  VectorJ X, u, p;            // 状态、输入、动量
  KalmanFilter *filter;       // ★ 唯一用裸 new 的地方
```

#### 构造函数逐行

```cpp
DKalmanObserver::DKalmanObserver(RobotDynamics* rd, MatrixJ& s, MatrixJ& h, MatrixJ& q, MatrixJ& r)
  : ExternalObserver(rd, ID_DKalmanObserver)
  , H(h)
  , X(VectorJ(2*jointNo))      // ← 状态是 2n 维
  , u(VectorJ(jointNo))
  , p(VectorJ(jointNo))
  , filter(0)
{
  MatrixJ A(2*jointNo, 2*jointNo), B(2*jointNo, jointNo), C(jointNo, 2*jointNo);
  B.setZero();
  C.setZero();
  for(int i = 0; i < jointNo; i++) {
    B(i, i) = 1;               // B = [I; 0]
    C(i, i) = 1;               // C = [I  0]
  }
  A.setZero();
  A.block(      0, jointNo, jointNo, jointNo) = h;   // A 右上 = H
  A.block(jointNo, jointNo, jointNo, jointNo) = s;   // A 右下 = S
  //  A 左上、左下保持 0 ✓（有 setZero）
  filter = new KalmanFilter(A, B, C);
  filter->setCovariance(q, r);
}
```

**逐行对照矩阵**：

$$
A = \begin{bmatrix} \mathbf{0} & H \\ \mathbf{0} & S\end{bmatrix}
\;\longleftrightarrow\;
\texttt{A.setZero(); A.block(0,n,n,n)=h; A.block(n,n,n,n)=s;}
$$

$$
B = \begin{bmatrix} I \\ \mathbf 0\end{bmatrix}
\;\longleftrightarrow\;
\texttt{B.setZero(); for(i<n) B(i,i)=1;}
$$

$$
C = \begin{bmatrix} I & \mathbf 0\end{bmatrix}
\;\longleftrightarrow\;
\texttt{C.setZero(); for(i<n) C(i,i)=1;}
$$

> **注意 `A.setZero()` 是有的** —— 对比 §10.1 里 `KalmanFilterContinous` 的 `AQ` 忘了 `setZero()`，这里做对了。

`filter = new KalmanFilter(A, B, C);` —— 全库**唯一用裸 `new` 的地方**，因此这三个 Kalman 观测器是仅有的需要写析构函数的类：

```cpp
DKalmanObserver::~DKalmanObserver() { delete filter; }
```

> 没有禁用拷贝构造/赋值 —— 拷贝一个 `DKalmanObserver` 会导致 **double-free**。实践中不会这么用（都是 `new` 出来存指针），但这是个潜在陷阱。改成 `std::unique_ptr<KalmanFilter>` 就自动安全了。

#### 每步逐行

```cpp
VectorJ DKalmanObserver::getExternalTorque(VectorJ& q, VectorJ& qd, VectorJ& tau, double dt)
{
  p = dyn->getM(q) * qd;                          // ① y = p = M·q̇（"测量"）
  u = tau - dyn->getG(q) - dyn->getFriction(qd);  // ② u = τ - G - τ_f
  u += dyn->getC(q, qd).transpose() * qd;         // ③ u += Cᵀq̇   → u = τ - β

  if(isRun) {
    X = filter->step(u, p, dt);                   // ④ ★ 3 参数版 → A 按连续处理
  } else {
    X.setZero();
    for(int i = 0; i < jointNo; i++) X(i) = p(i); // ⑤ X₀ = [p(0); 0]
    filter->reset(X);                             // ⑥ 同时 P₀ = 0
    isRun = true;
  }
  p = H * X.block(jointNo, 0, jointNo, 1);        // ⑦ ★ τ_ext = H·ω，复用 p
  return p;
}
```

**关键行解释**：

**②③ `u` 的构造**：$u = \tau - \tau_f - G + C^T\dot q = \tau - \tau_f - \beta$ —— 正是核心方程里的 $\tau - \beta$（扣掉摩擦）。分两行写是为了避免一个超长表达式。

**④ 用 3 参数版 `step()`**：传进 `KalmanFilter` 的 `A` 是**连续时间**矩阵，配 `step(u,y,dt)`（内部 `At = I + dt*A`）—— **自洽 ✓**。参见 §5.4 的坑。

**⑤ 首次不调用 `step()`**：第一帧的测量只用于初始化 $X_0$，不做滤波更新。所以首次返回 $H\cdot 0 = 0$ ✓。

**⑦ `X.block(jointNo, 0, jointNo, 1)`**：取状态向量的后 $n$ 个元素，即 $\omega$。Eigen 的 `block(startRow, startCol, rows, cols)` 对向量来说 `startCol=0, cols=1`。等价的更清晰写法是 `X.tail(jointNo)`。

#### 首次调用路径追踪

```
① p = M·q̇ = p(0)
② u = τ - G - τ_f
③ u += Cᵀq̇
⑤ X.setZero(); X(0..n-1) = p    → X₀ = [p(0); 0]
⑥ filter->reset(X)              → X = X₀, P = 0
⑦ p = H · X(n..2n-1) = H·0 = 0  ★
   return 0  ✓
```

初始状态 $X_0 = [p(0);\,0]$ 的含义：动量部分用实测值（准），扰动部分猜 0（未必准，但 $P_0=0$ 后每步预测都 `+Q`，不确定性会迅速注入，见 §5.4）。

#### 卡尔曼递推在这个系统上的具体形式

`filter->step(u, p, dt)` 里发生了什么（§5.4 是通用形式，这里代入本系统的 $A,B,C$ 看具体长相）。离散化 $A_d = I + \Delta t A = \begin{bmatrix}I & \Delta t H\\ 0 & I + \Delta t S\end{bmatrix}$，协方差按 $n\times n$ 分块 $P = \begin{bmatrix}P^{pp} & P^{p\omega}\\ P^{\omega p} & P^{\omega\omega}\end{bmatrix}$。

**预测步**（前向欧拉）：

$$
\begin{aligned}
\hat p^- &= \hat p + \Delta t\,(H\hat\omega + u) \\
\hat\omega^- &= (I + \Delta t\,S)\,\hat\omega \\
P^- &= A_d\,P\,A_d^T + Q
\end{aligned}
$$

第一式就是 $\dot p = u + H\omega$ 的一步积分：**用当前的扰动估计 $\hat\omega$ 去外推动量**。

**更新步**。因为 $C = [I\ 0]$，所以 $C\hat X^- = \hat p^-$，新息（innovation）是：

$$
\nu = y - C\hat X^- = \underbrace{p^{meas}}_{M(q)\dot q} - \hat p^-
$$

新息协方差和卡尔曼增益（注意 $CP^-C^T$ 就是取 $P^-$ 的左上块 $P^{-,pp}$）：

$$
S_{\nu} = P^{-,pp} + R, \qquad
K = P^- C^T S_\nu^{-1} = \begin{bmatrix}P^{-,pp}\\ P^{-,\omega p}\end{bmatrix} S_\nu^{-1}
= \begin{bmatrix}K_p\\ K_\omega\end{bmatrix}
$$

状态修正：

$$
\begin{bmatrix}\hat p\\ \hat\omega\end{bmatrix}
= \begin{bmatrix}\hat p^-\\ \hat\omega^-\end{bmatrix}
+ \begin{bmatrix}K_p\\ K_\omega\end{bmatrix}\nu
$$

##### 新息的物理意义（这个方法的核心直觉）

$$
\nu = p^{meas} - \hat p^- = \text{「实测动量」} - \text{「模型（含当前扰动估计）预测的动量」}
$$

就是**动量残差** —— 机器人实际的动量比模型预期多了多少。这份"多出来的动量"是未建模力矩的证据，被增益的两块分别处理：

- $K_p\nu$ 修正动量估计 $\hat p$
- $K_\omega\nu$ 修正**扰动估计** $\hat\omega$ ← 这才是我们要的外力

**关键**：$K_\omega = P^{-,\omega p}S_\nu^{-1}$ 里的**互协方差 $P^{\omega p}$** 是把动量残差引到扰动估计的桥梁。若 $P^{\omega p}=0$，扰动永远不会被更新。这个互协方差之所以会长出来，正是因为 $A_d$ 的非对角块 $\Delta t H$ —— 它在预测步把 $\omega$ 的不确定性泵进 $p$，从而在 $p$ 的残差和 $\omega$ 之间建立统计关联。**$H$ 既定义了扰动如何影响力矩，也决定了扰动能否被观测到。**

##### 可观性：什么时候扰动能被重构

$\omega$ 藏在状态里、只能通过 $p$ 间接观测，那它一定能被重构吗？不一定。看 $(A,C)$ 的可观性矩阵：

$$
C = [I\ 0], \quad
CA = [0\ H], \quad
CA^2 = [0\ HS], \ \dots
\;\Longrightarrow\;
\mathcal{O} = \begin{bmatrix} I & 0\\ 0 & H\\ 0 & HS\\ \vdots\end{bmatrix}
$$

上半部分 $[I\ 0]$ 已经张成动量子空间；能否张满整个 $2n$ 维，取决于右下角 $\begin{bmatrix}H\\ HS\\ HS^2\\ \vdots\end{bmatrix}$ 是否满秩 —— **这恰好是 $(S, H)$ 的可观性条件**。

$$
\boxed{\;\text{增广系统可观} \iff (S,\,H) \text{ 可观}\;}
$$

推论：
- **$H = I$（满秩）→ 无条件可观**，与 $S$ 无关。这就是为什么 $H=I$ 是安全默认值。
- **$H$ 奇异**（比如只想观测部分关节的外力）→ 必须验证 $(S,H)$ 可观，否则某些扰动分量估计不出来。

##### 与动量观测器的关系

取 $S=0,\ H=I$（默认配置），外部模型退化成 $\dot\omega = 0$、$\tau_{ext} = \omega$。此时扰动更新 $\hat\omega \leftarrow \hat\omega^- + K_\omega\nu$ 在结构上和动量观测器的

$$
\dot r = K_O(\tau_{ext} - r)
$$

**是一回事** —— 都是「用动量残差驱动外力估计朝真值收敛」。差别只在增益：

| | 增益来源 |
|---|---|
| `MomentumObserver` | 手调的**固定** $K_O$ |
| `DKalmanObserver` | 卡尔曼增益 $K_\omega = P^{\omega p}S_\nu^{-1}$，**由噪声统计自动整定、随 $P$ 时变**，稳态收敛到 Riccati 方程的最优解 |

所以可以把这个观测器理解成**「会自己调增益的动量观测器」** —— 代价是你得给出 $Q, R$（噪声统计）而不是直接给增益。

##### 无偏性交叉验证

设外力恒定 $\tau_{ext} = d^{*}$，用默认 $S=0, H=I$ 建模。则真实状态 $[p;\,\omega] = [p;\,d^{*}]$ 与模型自洽（$\omega \equiv d^{*}$ 满足 $\dot\omega = 0$）。因系统可观（$H=I$）且模型正确，卡尔曼估计渐近无偏：

$$
\hat\omega \longrightarrow d^{*} \quad\Longrightarrow\quad \hat\tau_{ext} = H\hat\omega \longrightarrow d^{*}
$$

**对常值外力无稳态偏差 ✓**。反之若把 $S$ 设错（比如给了衰减动态 $S=-\lambda I$ 却去测常值外力），模型和真实不符，会引入稳态偏差 —— 这是选 $S$ 时要当心的。

#### 调参

```matlab
S = zeros(JNT,JNT);                              % 扰动动态
H = eye(JNT,JNT);                                % 扰动 → 力矩映射
Q = blkdiag(0.002*eye(JNT), 0.3*eye(JNT));      % 过程噪声 (2n×2n)
R = 0.05*eye(JNT,JNT);                           % 测量噪声 (n×n)
```

**$S = 0$ → 随机游走模型**：

$$
\dot\omega = 0 \cdot \omega = 0 \quad\text{（+ 过程噪声）}
$$

即假设外力是**分段常值**，变化完全由过程噪声驱动。这是最常用的选择（测试代码全用它）。如果你知道外力的动态特性（比如已知频率 $\omega_0$ 的振动），可以设 $S = \begin{bmatrix}0 & 1\\ -\omega_0^2 & 0\end{bmatrix}$ 这类结构。

**$H = I$** → 扰动状态直接就是外部力矩。

**$Q$ 的分块结构是主要旋钮**：

$$
Q = \begin{bmatrix} 0.002\,I & 0 \\ 0 & 0.3\,I\end{bmatrix}
$$

- 前 $n$ 个对应动量 $p$ —— **小**（0.002），因为动力学模型比较可信
- 后 $n$ 个对应扰动 $\omega$ —— **大**（0.3），因为外力可以任意突变

**这个比例（$0.3/0.002 = 150$）决定了响应速度 vs 噪声抑制**。比例越大越信任「扰动会突变」这个假设，响应越快但越容易被噪声骗。

**$R$** 是动量"测量"的噪声。$R$ 越大越不信测量，输出越平滑但越滞后。

> 卡尔曼滤波的调参本质上只看 $Q/R$ 的**比值**，同时放大两者结果不变。

> ⚠️ **$Q$ 是「每步」协方差，和 $\Delta t$ 绑定**。3 参数版 `step(u,y,dt)` 里 $A,B$ 都按 $\Delta t$ 缩放了，唯独 `P = At*P*At' + Q` 的 $Q$ 没有（§10.4）。所以**控制周期一改，这些 $Q$ 值必须重新整定**。`DKalmanObserverExp` 用 Van Loan 正确离散化了 $Q_d$，所以它的 $Q$ 才是「每秒」的物理量、与周期无关 —— 这也是两者参数量级差 100 倍的原因之一（§6.5）。

#### RNEA 版差异

[disturbance_kalman_filter_rnea.h:104-106](lib/disturbance_kalman_filter_rnea.h#L104-L106)：

```cpp
p = dyn->rnea(q, zero, qd);                                        // M·q̇
u = tau - dyn->rnea(q, zero, zero, GRAVITY) - dyn->getFriction(qd);// τ - G - τ_f
u += dyn->tranCqd(q, qd);                                          // + Cᵀq̇
```

---

### 6.5 [lib/disturbance_kalman_filter_exp.h](lib/disturbance_kalman_filter_exp.h) —— 卡尔曼观测器 · 精确离散版

和 `DKalmanObserver` **结构完全相同**，只有三处差异：

```cpp
// ① 多一个成员，缓存 M
MatrixJ H, M;

// ② 换成 Van Loan 精确离散化的滤波器
filter = new KalmanFilterContinous(A, B, C);

// ③ 每步更新位形相关的 R
M = dyn->getM(q);
p = M * qd;
...
filter->updateR(M);        // R ← M·R·Mᵀ
```

#### 逐行代码（只看差异部分）

```cpp
VectorJ DKalmanObserverExp::getExternalTorque(VectorJ& q, VectorJ& qd, VectorJ& tau, double dt)
{
  M = dyn->getM(q);                               // ① 缓存 M（原版直接用不存）
  p = M * qd;
  u = tau - dyn->getG(q) - dyn->getFriction(qd);
  u += dyn->getC(q, qd).transpose() * qd;

  filter->updateR(M);                             // ② ★ 无条件调用，isRun 之前

  if(isRun) {
    X = filter->step(u, p, dt);
  } else {
    X.setZero();
    for(int i = 0; i < jointNo; i++) X(i) = p(i);
    filter->reset(X);
    isRun = true;
  }
  p = H * X.block(jointNo, 0, jointNo, 1);
  return p;
}
```

**① 为什么要缓存 `M`**：原版里 `dyn->getM(q)` 只用一次（`p = dyn->getM(q) * qd`），这里要用两次（算 `p` 和传给 `updateR`），所以存起来避免重复计算。**这个细节体现了作者对性能的注意。**

**② `updateR(M)` 在 `if(isRun)` 之前**：即使首次调用也会更新 $R$。合理 —— `reset()` 后的第一次 `step()` 就需要正确的 $R$。

#### 理论改进（见 §5.5）


|            | `DKalmanObserver`        | `DKalmanObserverExp`                          |
| ---------- | ------------------------ | --------------------------------------------- |
| 离散化     | 前向欧拉$I + \Delta t A$ | Van Loan 矩阵指数                             |
| $Q$ 离散化 | 无（原样用）             | $\int_0^{\Delta t}e^{A\tau}Qe^{A^T\tau}d\tau$ |
| $R$ 离散化 | 无                       | $R_d = \dfrac{MRM^T}{\Delta t}$               |
| $R$ 随位形 | ❌                       | ✅                                            |

#### 参数量级完全不同

```cpp
// DKalmanObserver
Q(0,0) = 0.002;  Q(1,1) = 0.002;  Q(2,2) = 0.3;  Q(3,3) = 0.3;    R *= 0.05;
// DKalmanObserverExp
Q(0,0) = 0.2;    Q(1,1) = 0.2;    Q(2,2) = 30;   Q(3,3) = 30;     R *= 0.0005;
```

$Q$ 大 100 倍，$R$ 小 100 倍 —— 因为离散化方式变了（$Q_d$ 现在会被 $\Delta t$ 积分缩小，$R_d$ 会被 $\Delta t$ 除大）。

> **参数绝对不能在两者之间照搬。**

> ⚠️ 这个观测器依赖 `KalmanFilterContinous`，而后者有**未初始化内存**的 bug（§10.1）。**在修掉之前，这个观测器的结果不可信。** 这或许能解释为什么 [matlab/test.m:106-126](matlab/test.m#L106-L126) 里它是**唯一被整段注释掉**的观测器。

---

### 6.6 [lib/filtered_dyn_observer.h](lib/filtered_dyn_observer.h) —— 滤波动力学观测器

**最优雅的一个**，也是唯一完全没有反馈的。

#### 理论 —— 用滤波器代数吸收微分

从核心方程直接解出：

$$
\tau_{ext} = \dot p + \big(\tau_f + G - C^T\dot q - \tau\big)
$$

$\dot p$ 算不了。**但**：碰撞检测本来就要对结果滤波，那不如**只求一个低通滤波后的估计**：

$$
\hat\tau_{ext} = W(s)\,\tau_{ext}, \qquad W(s) = \frac{\omega}{s+\omega}
$$

展开：

$$
\hat\tau_{ext} = \underbrace{W(s)\,\dot p}_{\text{第一项}} + W(s)\big(\tau_f + G - C^T\dot q - \tau\big)
$$

第一项在 $s$ 域是 $W(s)\cdot s\cdot p$，里面的 $s$ 就是微分。**关键的代数恒等式**：

$$
\boxed{\;\frac{\omega\,s}{s+\omega} = \frac{\omega(s + \omega - \omega)}{s+\omega} = \omega - \frac{\omega^2}{s+\omega}\;}
$$

于是：

$$
W(s)\,\dot p = \underbrace{\omega\,p}_{\text{纯代数}} \;+\; \underbrace{\left(\frac{-\omega^2}{s+\omega}\right)p}_{\text{稳定的一阶滤波器}}
$$

**微分被消掉了！** 右边只剩下 $p$ 的代数放大和一个**稳定、真分式**的滤波器。

$$
\boxed{\;\hat\tau_{ext} = \omega\,p + F_2(s)\,p + F_1(s)\big(\tau_f + G - C^T\dot q - \tau\big)\;}
$$

其中 $F_1(s) = \dfrac{\omega}{s+\omega}$，$F_2(s) = \dfrac{-\omega^2}{s+\omega}$。

> **这就是 `FilterF2` 存在的全部理由** —— 它那个奇怪的传递函数 $-\omega^2/(s+\omega)$ 不是凭空来的，是从这个恒等式推出来的。§5.3 里看不懂它为什么长这样，到这里就通了。

#### 成员变量表

```cpp
private:
  FilterF1 f1;      // ω/(s+ω)
  FilterF2 f2;      // -ω²/(s+ω)
  VectorJ p, res;
```

**只有 2 个滤波器 + 2 个向量** —— 全库最简单的观测器。没有任何反馈状态。

#### 逐行代码

```cpp
VectorJ FDynObserver::getExternalTorque(VectorJ& q, VectorJ& qd, VectorJ& tau, double dt)
{
  p = dyn->getM(q) * qd;                              // ① p = M·q̇

  if(isRun) {
    res = f2.filt(p, dt) + f2.getOmega() * p;         // ② W(s)·ṗ = F₂(s)p + ω·p
    //    └── F₂(s)·p ──┘   └──── ω·p ────┘
    p = dyn->getFriction(qd) + dyn->getG(q) - dyn->getC(q, qd).transpose() * qd;
    //  └──────────── τ_f + G - Cᵀq̇ ────────────┘                            ③
    p -= tau;                                         // ④ τ_f + G - Cᵀq̇ - τ
    res += f1.filt(p, dt);                            // ⑤ + F₁(s)·(...)
  } else {
    f2.set(p);                                        // ⑥ 稳态初始化 F₂
    p = dyn->getFriction(qd) + dyn->getG(q) - dyn->getC(q, qd).transpose() * qd;
    p -= tau;
    f1.set(p);                                        // ⑦ 稳态初始化 F₁
    res.setZero();                                    // ⑧ 首次输出 0
    isRun = true;
  }
  return res;
}
```

**关键行解释**：

**② 对照公式**：`f2.getOmega()` 返回 $\omega$（截止频率），所以这一行就是

$$
\texttt{res} = F_2(s)\,p + \omega\,p = W(s)\,\dot p
$$

§5.3 提到 `getOmega()` 命名容易误导（返回 `cut` 而不是成员 `omega`）—— 现在能看出它返回的确实是推导里需要的 $\omega$。

**③④ `p` 被复用**：从动量变成 $\tau_f + G - C^T\dot q - \tau$。注释写着 `// reuse`。

**⑥⑦ 稳态初始化的作用**：

首次调用时 `f2.set(p)` 让 $F_2$ 的输出立刻等于稳态值 $-\omega p$，于是第二次调用时：

$$
\texttt{res} \approx -\omega p + \omega p = 0
$$

**没有启动瞬态** —— 如果不做这个初始化，滤波器从零状态启动会产生一个巨大的虚假尖峰（因为 $p$ 从"滤波器认为的 0"跳到实际值），足以触发误报警。

同理 `f1.set(p)` 让 $F_1$ 的输出立刻等于其输入的稳态值。

**⑧ 首次显式返回 0**：即使有 ⑥⑦ 的初始化，第一帧仍显式清零，保证 §4.3 的不变量。

#### 首次调用路径追踪

```
① p = M·q̇
⑥ f2.set(p)            → F₂ 状态 = 稳态（y = -ω·p）
   p = τ_f + G - Cᵀq̇ - τ
⑦ f1.set(p)            → F₁ 状态 = 稳态（y = p）
⑧ res.setZero()        ★
   return 0  ✓
```

#### 特点与取舍


|    |                                                               |
| -- | ------------------------------------------------------------- |
| ✅ | **只有一个参数**（截止频率 $\omega$），最好调                 |
| ✅ | **无状态反馈** → 结构上不可能发散                            |
| ✅ | 便宜：1 次 M/C/G，无矩阵求逆                                  |
| ✅ | 无启动瞬态（稳态初始化做得很到位）                            |
| ❌ | **纯前馈** → 对模型误差零抑制能力，误差 $1{:}1$ 直接进入输出 |
| ❌ | 精度完全取决于动力学参数的准确性                              |

**这正是它和 `MomentumObserver` 的本质区别**：

$$
\text{MO:}\quad \dot r = K_O(\tau_{ext} - r) \qquad \text{有反馈回路，能一定程度抑制模型误差}
$$

$$
\text{FDyn:}\quad \hat\tau_{ext} = W(s)\,\tau_{ext} \qquad \text{纯开环，模型错多少输出就错多少}
$$

> **这也解释了为什么需要 `dynamic_calibration`** —— 参数不准的话，FDyn 直接废掉。

#### 调参

```cpp
FDynObserver fd_observer(&robot, 8, TSTEP);   // ω = 8 rad/s, T = 0.01 s
```

只有截止频率 $\omega$ 一个参数（$\Delta t$ 由调用方每步传入，构造时那个只是初值）：

- $\omega$ 大 → 响应快，但噪声和模型误差都进来了
- $\omega$ 小 → 平滑，但碰撞检测延迟大

$\omega = 8$ rad/s $\approx 1.27$ Hz，对应时间常数 125 ms。

> ⚠️ **`8` 的单位是 rad/s，不是 Hz**，尽管参数名叫 `cutOffHz`；`0.01` 是采样**周期（秒）**，尽管参数名叫 `sampHz`。见 §10.3。

#### RNEA 版差异

[filtered_dyn_observer_rnea.h:44,48](lib/filtered_dyn_observer_rnea.h#L44)：

```cpp
p = dyn->rnea(q, zero, qd);                                                    // M·q̇
...
p = dyn->getFriction(qd) + dyn->rnea(q, zero, zero, GRAVITY) - dyn->tranCqd(q, qd);
//                         └────── G ──────┘                   └── Cᵀq̇ ──┘
```

---

### 6.7 [lib/filtered_range_observer.h](lib/filtered_range_observer.h) —— 范围观测器

**最容易被误解的一个。它不估计外部力矩 —— 它估计阈值。**

#### 和 FDynObserver 的两处差异

对比两者代码，逻辑几乎逐字相同，但：

```cpp
VectorJ FRangeObserver::getExternalTorque(VectorJ& q, VectorJ& qd, VectorJ& tau, double dt)
{
  p = dyn->getM(q) * qd;
  p *= shift;                       // ① ★ 所有模型项乘以 k

  if(isRun) {
    res = f2.filt(p, dt) + f2.getOmega() * p;
    p = dyn->getFriction(qd) + dyn->getG(q) - dyn->getC(q, qd).transpose() * qd;
    p *= shift;                     // ① ★ 同上
    res += f1.filt(p, dt);          // ② ★ 完全没有 tau！
  } else {
    f2.set(p);
    p = dyn->getFriction(qd) + dyn->getG(q) - dyn->getC(q, qd).transpose() * qd;
    p *= shift;
    f1.set(p);
    res.setZero();
    isRun = true;
  }
  return res;
}
```

**`FDynObserver` 里有 `p -= tau;`，这里没有。这不是遗漏，是设计。**

#### 它在算什么

`FRangeObserver` 的输出是：

$$
\text{res} = k \cdot \Big[\omega\,p + F_2(s)\,p + F_1(s)\big(\tau_f + G - C^T\dot q\big)\Big]
$$

对比 `FDynObserver`：

$$
\hat\tau_{ext} = \omega\,p + F_2(s)\,p + F_1(s)\big(\tau_f + G - C^T\dot q - \tau\big)
$$

**差异**：乘了 $k$，且**不含 $\tau$**。

**物理含义**：假设动力学模型有 $\pm k$ 的**相对误差**，$\hat\tau_{ext}$ 会漂多少。

- **乘 $k$**：模型项按比例缩放 —— 若 $M, C, G, \tau_f$ 全都有 $k$ 的相对误差，输出的误差就是模型项的 $k$ 倍
- **不含 $\tau$**：因为 $\tau$ 是**测量值**，不受动力学参数误差影响 —— 模型参数错了不会让力矩读数变化

#### 用法：动态阈值

[matlab/test.m:147-169](matlab/test.m#L147-L169) 揭示了真实用途：

```matlab
delta = 0.1;
id_dr1 = calllib('observers','configFilterRangeObserver',-1,cutOff,timeStep, delta);  % +10%
id_dr2 = calllib('observers','configFilterRangeObserver',-1,cutOff,timeStep,-delta);  % -10%

for i = 2:size(cur,1)
    calllib('observers','getExternalTorque',id_dr1,res,q(i,:),qd(i,:),K.*cur(i,:),tm(i)-tm(i-1));
    ext_up(i,:) = res.Value(:);
    calllib('observers','getExternalTorque',id_dr2,res,q(i,:),qd(i,:),K.*cur(i,:),tm(i)-tm(i-1));
    ext_low(i,:) = res.Value(:);
end

n = 1;  % joint index
plot(tm,[ext_up(:,n), ext_low(:,n), ext_df(:,n)]);
title("Range");
```

开**两个实例**（$+k$ 和 $-k$），画出 `FDynObserver` 输出的**上下包络带**：

```
      τ̂_ext (FDyn)
         │      ╱╲          ← 真实碰撞：冲出包络带 → 报警
    ─────┼─────╱──╲──────   ← ext_up   (+10% 模型误差带)
         │   ╱      ╲
    ─────┼──────────────    ← ext_low  (-10% 模型误差带)
         └──────────────→ t
```

判据：

$$
\hat\tau_{ext} \notin \big[\text{ext\_low},\ \text{ext\_up}\big] \quad\Longrightarrow\quad \text{碰撞}
$$

**这就是「自适应阈值」碰撞检测**：不用固定阈值（高速运动时模型误差大、容易误报；静止时阈值又过于保守），而是用一个**随位形和速度变化**的置信带。

> 这是论文标题里 **"Practical Aspects"** 的精髓 —— 理论上七个观测器都能估 $\tau_{ext}$，但**工程上你还需要知道什么时候该相信它**。

#### 调参

```cpp
FRangeObserver fr_observer(&robot, 8, TSTEP, -0.1);
//                                  ω   T    k = -10%
```

$k$ 就是你对动力学参数相对误差的估计。做过 `dynamic_calibration` 且验证 `rre` 很低的话可以取小（5%）；参数是从 URDF 抄的话得取大（20%+）。

$\omega$ 和 $T$ **必须和配套的 `FDynObserver` 完全一致**，否则包络带和被包络的信号不在同一个频带上，比较没有意义。

> **没有 RNEA 版本。** 要在 RNEA 路线上用自适应阈值，得自己照着 `filtered_dyn_observer_rnea.h` 写一个。

---

## 7. 观测器选型对比

### 7.1 计算成本


| 观测器                | 动力学调用       | 矩阵求逆               | 其他                                  | 相对成本      |
| --------------------- | ---------------- | ---------------------- | ------------------------------------- | ------------- |
| `MomentumObserver`    | 1×M, 1×C, 1×G | 无                     | 梯形积分                              | ★ 最低       |
| `FDynObserver`        | 1×M, 1×C, 1×G | 无                     | 2 个 IIR + 2 次`tan()`                | ★ 低         |
| `FRangeObserver`      | 1×M, 1×C, 1×G | 无                     | 同上                                  | ★ 低         |
| `SlidingModeObserver` | 1×M, 1×C, 1×G | 无                     | $n$ 次 `tanh()` + $n$ 次 `sqrt()`     | ★★ 中       |
| `DKalmanObserver`     | 1×M, 1×C, 1×G | 1× ($n{\times}n$)     | $2n{\times}2n$ 矩阵乘若干             | ★★★ 高     |
| `DisturbanceObserver` | 1×M, 1×C, 1×G | **2× ($n{\times}n$)** | —                                    | ★★★ 高     |
| `DKalmanObserverExp`  | 1×M, 1×C, 1×G | 1× ($n{\times}n$)     | **2× 矩阵指数**（各 5 次矩阵乘）+ KF | ★★★★ 最高 |

RNEA 版本一律再乘约 10 倍（§5.2）。

### 7.2 特性对比


| 观测器                | 参数个数         | 收敛性         | 抗模型误差       | 数值稳定性        | 积分方法         |
| --------------------- | ---------------- | -------------- | ---------------- | ----------------- | ---------------- |
| `MomentumObserver`    | $n$ ($k$)        | 渐近（一阶）   | 中（有反馈）     | 好                | 梯形             |
| `DisturbanceObserver` | 3→1 ($k$)       | 渐近（一阶）   | 中               | **很好**          | 后向欧拉（隐式） |
| `SlidingModeObserver` | $4n \to 2n$      | **有限时间**   | 好               | 一般（可能抖振）  | 前向欧拉（显式） |
| `DKalmanObserver`     | $Q, R$           | 随机最优       | 好（可建模）     | 好                | 前向欧拉         |
| `DKalmanObserverExp`  | $Q, R$           | 随机最优       | 好               | 好（但见 §10.1） | Van Loan         |
| `FDynObserver`        | **1** ($\omega$) | 即时（无动态） | **差**（无反馈） | **最好**          | 无（纯滤波）     |
| `FRangeObserver`      | 2 ($\omega, k$)  | —             | —               | 最好              | 无               |

### 7.3 该用哪个

- **先用 `MomentumObserver`** —— 最经典，最便宜，一个参数量级好猜（30~50），出问题最好排查。**绝大多数场合够用。**
- **`FDynObserver`** —— 如果你的动力学参数很准（跑过 `dynamic_calibration` 且 `rre` 很低），它最省心：一个参数，不会发散，无启动瞬态。
- **`FDynObserver` + `FRangeObserver`×2** —— 想做自适应阈值碰撞检测，这是**论文推荐的组合**。
- **`DKalmanObserver`** —— 如果你对噪声统计有把握，或想把外力的先验动态（$S$ 矩阵）建进去。
- **`SlidingModeObserver`** —— 需要有限时间收敛的保证（安全认证场合）。准备好调参和处理抖振。
- **`DisturbanceObserver`** —— $\Delta t$ 抖动大、其他观测器不稳定时，隐式离散救场。
- **`DKalmanObserverExp`** —— **先修 §10.1 的 bug**。

---

## 8. 测试与 MATLAB 封装

### 8.1 [tests/](tests/) —— 是示例，不是测试

**这里没有断言，没有自动判定。** 全部是「跑一遍 → 存 CSV → gnuplot 画图 → 人眼看」。README 也只说 "Tests / examples"。

```bash
make test     # M/C/G 版观测器 → force.csv + input.csv → gnuplot
make rnea     # RNEA 版观测器  → force.csv → gnuplot
make kalman   # 卡尔曼滤波器单独验证 → estimation.csv → gnuplot
make dyn      # 编译成 libobservers.so
```

对比 `dynamic_calibration` 的 `tests/`（那里有真 `assert`，还拿 MATLAB Robotics Toolbox 做独立参照）—— 这个库的测试薄弱得多。

#### 玩具模型：[tests/double_link.h](tests/double_link.h)

标准的平面 2 连杆机械臂，参数写死：

```cpp
m1 = 1;     m2 = 1;       // 质量 (kg)
l1 = 0.5;   l2 = 0.5;     // 连杆长度 (m)
lc1 = 0.25; lc2 = 0.25;   // 质心位置 (m)
I1 = 0.3;   I2 = 0.2;     // 转动惯量 (kg·m²)
```

惯性矩阵（教科书标准形式）：

$$
M(q) = \begin{bmatrix}
m_1l_{c1}^2 + m_2(l_1^2 + l_{c2}^2 + 2l_1l_{c2}\cos q_2) + I_1 + I_2 & m_2(l_{c2}^2 + l_1l_{c2}\cos q_2) + I_2 \\
m_2(l_{c2}^2 + l_1l_{c2}\cos q_2) + I_2 & m_2l_{c2}^2 + I_2
\end{bmatrix}
$$

科氏矩阵（令 $h = -m_2l_1l_{c2}\sin q_2$）：

$$
C(q,\dot q) = \begin{bmatrix} h\dot q_2 & h(\dot q_1 + \dot q_2) \\ -h\dot q_1 & 0 \end{bmatrix}
$$

```cpp
MatrixJ DoubleLink::getC(VectorJ& q, VectorJ& qd)
{
  double h = -m2*l1*lc2*sin(q(1));
  C(0,0) = h*qd(1);
  C(0,1) = h*(qd(0)+qd(1));
  C(1,0) = -h*(qd(0));
  return C;                       // ← C(1,1) 从未赋值！
}
```

> **`C(1,1)` 从未赋值** —— 靠构造函数的 `C.setZero()` 保持为 0。这是**对的**（标准 2 连杆模型 $C_{22} = 0$），但依赖构造时的初始化，读的时候容易以为漏了。**如果有人后来加了个 `C.setIdentity()` 之类的，这里会静默出错。**

重力项：

$$
G(q) = \begin{bmatrix}(m_1l_{c1} + m_2l_1)g\cos q_1 + m_2l_{c2}g\cos(q_1+q_2) \\ m_2l_{c2}g\cos(q_1+q_2)\end{bmatrix}
$$

**摩擦被完全忽略**：

```cpp
VectorJ getFriction(VectorJ& qd) { return fric; }   // fric 在构造后就是零向量，再没动过
```

> 真机上摩擦是**最大的误差源**（`dynamic_calibration` 的 README 花了大篇幅讲摩擦辨识）。所以这些测试结果**比真实情况乐观得多**。

#### 测试怎么造外力：[tests/test.cpp:89-95](tests/test.cpp#L89-L95)

```cpp
tau = robot.getM(q)*q2d + robot.getC(q,qd)*qd + robot.getG(q);   // 理想力矩
#ifdef SET_TORQUE
    if(t > 1 && t < 2) {      // 1~2 秒之间
      tau(0) -= 0.5;          // ★ 人为扣掉 0.5 N·m
      tau(1) -= 0.5;
    }
#endif
```

**符号逻辑要想清楚**：

1. 正向动力学算出**无外力时**该有的 $\tau$
2. 然后**减去** 0.5

观测器看到的是「力矩比模型预期少了 0.5」。根据动力学方程

$$
M\ddot q + C\dot q + G = \underbrace{\tau_{\text{给的}}}_{\tau_{\text{理想}} - 0.5} + \tau_{ext}
$$

$$
\Longrightarrow \tau_{ext} = \underbrace{M\ddot q + C\dot q + G}_{\tau_{\text{理想}}} - (\tau_{\text{理想}} - 0.5) = +0.5
$$

**正确的输出应该是 1~2 秒之间出现一个 $+0.5$ 的方波**，其余时间为 0。

轨迹是两个不同频率的正弦：

```cpp
#define OMEGA1 1.3
#define OMEGA2 0.8
q(0) = sin(1.3t);   qd(0) = 1.3cos(1.3t);   q2d(0) = -1.69sin(1.3t);
q(1) = sin(0.8t);   qd(1) = 0.8cos(0.8t);   q2d(1) = -0.64sin(0.8t);
```

不同频率保证两个关节不同步、持续激励（和 `dynamic_calibration` 的傅里叶激励轨迹是同一个思路）。

> **这是理想仿真**：$q, \dot q, \ddot q$ 解析给出，**无噪声**，模型和"真实"**完全一致**。所以它只能验证**观测器逻辑对不对**，完全不能反映真机性能。真机上噪声和模型误差才是主角。

#### `#ifdef` 一次只能开一个

```cpp
//#include "../lib/momentum_observer.h"
//#include "../lib/disturbance_observer.h"
#include "../lib/sliding_mode_observer.h"      // ← 当前只有这个启用
//#include "../lib/disturbance_kalman_filter.h"
```

如 §3.1 所说，这是 header-only 非 inline 实现的**结构性约束**，不是风格选择。要换观测器就改注释重新编译。

`test.cpp` 里所有观测器的配置和调用都包在 `#ifdef <头文件的 include guard>` 里：

```cpp
#ifdef SLIDING_MODE_OBSERVER_H      // ← 这个宏来自 sliding_mode_observer.h 的 include guard
  VectorJ T1(2), S1(2), T2(2), S2(2);
  ...
#endif
```

**很聪明的技巧** —— 用 include guard 当"这个头文件被包含了吗"的探测器，不需要额外的配置宏。

⚠️ 但那些被注释掉的分支里有**编译不过的死代码**，见 §10.2。

### 8.2 [matlab/](matlab/) —— C 接口封装

MATLAB 不能直接调 C++（名字修饰、类、模板），所以有一层 `extern "C"` 封装。

#### 设计：句柄表模式

```cpp
static DoubleLink robot;                             // 全局机器人对象
static ExternalObserver *observer[ARRAY_LEN] = {0};  // 观测器数组
static int _nextIndex = 0;                           // 下一个可用槽位
```

MATLAB 侧拿到的是**整型 ID**（数组下标），不是指针 —— 因为 MATLAB 处理不了 C++ 对象。

```cpp
#define ADD_NEW -1
// 传 -1 → 新建，返回新 ID；传已有 ID → 更新该观测器的参数
int configMomentumObserver(int index, double k[JOINT_NO]);
```

#### 每个 `config*` 的统一模式

```cpp
int configMomentumObserver(int ind, double *k)
{
  // ① 把 C 数组转成 Eigen 类型
  MomentumObserver* ptr;
  VectorJ vk(JOINT_NO);
  for(int i = 0; i < JOINT_NO; i++)
    vk(i) = k[i];

  // ② 新建路径
  if(ind == ADD_NEW) {
    if(_nextIndex == ARRAY_LEN) return ERR_NO_SLOTS;   // 槽位满
    ptr = new MomentumObserver(&robot, vk);
    observer[_nextIndex] = ptr;
    return _nextIndex++;                               // ★ 返回旧值再自增
  } else if(ind < ADD_NEW || ind >= ARRAY_LEN) {
    return ERR_WRONG_INDEX;                            // ③ 边界检查
  }

  // ④ 更新路径
  ptr = (MomentumObserver*) observer[ind];
  if(ptr->type() != ID_MomentumObserver)
    return ERR_WRONG_TYPE;                             // ★ 运行时类型检查
  ptr->settings(vk);
  return ind;
}
```

**关键点**：

- **② `return _nextIndex++`**：后缀自增，返回的是**自增前**的值（即刚占用的槽位），同时 `_nextIndex` 指向下一个空位。经典写法。
- **④ `type()` 的运行时检查**：防止把 ID 传错给另一种观测器的 config 函数。因为 C 接口丢失了类型信息，这个检查是必要的补偿。**这就是 §4.2 那个 ID 系统存在的理由。**
- **错误码**：`ERR_WRONG_INDEX = -1`、`ERR_NO_SLOTS = -2`、`ERR_WRONG_TYPE = -3`。注意 `ADD_NEW` 也是 `-1`，和 `ERR_WRONG_INDEX` **数值相同** —— 输入的 -1 表示"新建"，输出的 -1 表示"索引错误"。方向不同所以不冲突，但容易看晕。

#### `freeAll()` 的设计缺陷

```cpp
void freeAll()
{
  MomentumObserver *mo;  DisturbanceObserver *dis;  SlidingModeObserver *sm;
  DKalmanObserver *dk;   FDynObserver *fd;  FRangeObserver *fr;  DKalmanObserverExp *dke;
  for(int i = 0; i < _nextIndex; i++) {
    switch (observer[i]->type()) {
    case ID_MomentumObserver:    mo = (MomentumObserver*) observer[i];    delete mo;  break;
    case ID_DisturbanceObserver: dis = (DisturbanceObserver*) observer[i]; delete dis; break;
    case ID_SlidingModeObserver: sm = (SlidingModeObserver*) observer[i];  delete sm;  break;
    case ID_DKalmanObserver:     dk = (DKalmanObserver*) observer[i];      delete dk;  break;
    case ID_FDynObserver:        fd = (FDynObserver*) observer[i];         delete fd;  break;
    case ID_FRangeObserver:      fr = (FRangeObserver*) observer[i];       delete fr;  break;
    case ID_DKalmanObserverExp:  dke = (DKalmanObserverExp*) observer[i];  delete dke; break;
    }
    observer[i] = 0;
  }
  _nextIndex = 0;
}
```

**这整个 switch 是不必要的。** `ExternalObserverBase` 在 [external_observer.h:150](lib/external_observer.h#L150) 有虚析构：

```cpp
virtual ~ExternalObserverBase() = default;
```

所以

```cpp
for(int i = 0; i < _nextIndex; i++) { delete observer[i]; observer[i] = 0; }
_nextIndex = 0;
```

**一行就够**，虚析构会正确分派到派生类。

这个 switch 是 `7e4cf07 Add virtual destructors` 提交**之前**的遗留 —— 加了虚析构之后没有回来简化它。

**代价**：每加一个新观测器就得记得改这个 switch，忘了就**静默内存泄漏**（没有 `default` 分支，匹配不上就什么都不做，指针被 `observer[i] = 0` 抹掉，永久泄漏）。

#### 数据流：列主序转换

```cpp
int k = 0;
for(int c = 0; c < JOINT_NO; c++) {          // ← 列在外
    for(int r = 0; r < JOINT_NO; r++, k++) {  // ← 行在内
      mS(r,c) = S[k];
      mH(r,c) = H[k];
      mR(r,c) = R[k];
    }
}
```

**MATLAB 的矩阵是列主序展平的**：`A = [1 2; 3 4]` 传到 C 是 `{1, 3, 2, 4}`。所以外层循环列、内层循环行，才能正确还原。搞反了矩阵就转置了。

$Q$ 是 $2n\times 2n$，单独一个循环：

```cpp
k = 0;
for(int c = 0; c < 2*JOINT_NO; c++) {
    for(int r = 0; r < 2*JOINT_NO; r++, k++)
      mQ(r,c) = Q[k];
}
```

#### `getExternalTorque` 的封装

```cpp
int getExternalTorque(int ind, double* ext, double *q, double *qd, double *tau, double dt)
{
  if(ind < 0 || ind >= _nextIndex) return ERR_WRONG_INDEX;
  ExternalObserver* ob = observer[ind];

  VectorJ vext(JOINT_NO), vq(JOINT_NO), vqd(JOINT_NO), vtau(JOINT_NO);   // ← 每次都构造！
  for(int i = 0; i < JOINT_NO; i++) {
    vq(i) = q[i];  vqd(i) = qd[i];  vtau(i) = tau[i];
  }

  vext = ob->getExternalTorque(vq,vqd,vtau,dt);

  for(int i = 0; i < JOINT_NO; i++) ext[i] = vext(i);
  return 0;
}
```

> **讽刺的是**：库里所有观测器都费尽心思用成员变量避免堆分配（§5.1），而这个封装层**每次调用都新建 4 个 `VectorJ`**（各一次堆分配）。对 MATLAB 调用来说无所谓（MATLAB 的开销大得多），但如果你照抄这个模式做实时 ROS 节点，就把库的优化全浪费了。

#### [matlab/test.m](matlab/test.m) —— 论文的对比实验

这是全库最有价值的**用法示例**，七种观测器跑同一组数据：

```matlab
if not(libisloaded('observers'))
    loadlibrary('observers.so','observers.h')
end

data = csvread('input.csv');
tm   = data(:,1);      % 时间
q    = data(:,2:3);    % 角度
qd   = data(:,4:5);    % 速度
cur  = data(:,6:7);    % 电流/力矩

JNT = 2;
K = [1,1];    % ★ current to torque, N/A  ← dynamic_calibration 的产物
res = libpointer('doublePtr',zeros(1,JNT));   % 输出缓冲区

% ---- 每个观测器都是同一个模式 ----
id_mo = calllib('observers','configMomentumObserver',-1,Kmo);   % -1 = 新建
for i = 2:size(cur,1)
    calllib('observers','getExternalTorque', id_mo, res, ...
            q(i,:), qd(i,:), K.*cur(i,:), tm(i)-tm(i-1));
    %                        └─ 电流→力矩 ─┘  └── dt ──┘
    ext_mo(i,:) = res.Value(:);
end
...
calllib('observers','freeAll');   % 清内存
unloadlibrary('observers');
```

**关键行**：

- **`K.*cur(i,:)`** —— **电流乘以驱动增益变成力矩**。这一行就是两个项目的接合点。`test.m` 里 `K = [1,1]`（仿真数据本来就是力矩），真机上这里应该填 `estimate_drive_gains` 的输出。
- **`tm(i)-tm(i-1)`** —— 每步的实际 `dt`，从时间戳算。所以循环从 `i = 2` 开始。
- **`libpointer('doublePtr', ...)`** —— MATLAB 的输出参数机制。`res.Value` 取回内容。

`input.csv` 由 [tests/test.cpp](tests/test.cpp) 生成（它同时写 `force.csv` 和 `input.csv`）。

#### ⚠️ `matlab/` vs `tests/` 的 observers.cpp

**两个文件几乎同名，但不一样**：


|                               | `tests/observers.cpp` | `matlab/observers.cpp`         |
| ----------------------------- | --------------------- | ------------------------------ |
| 行数                          | 216                   | 303                            |
| `ARRAY_LEN`                   | 20                    | 30                             |
| 观测器数                      | 5                     | **7**（多 FRange、DKalmanExp） |
| `freeAll()` 重置 `_nextIndex` | ❌                    | ✅                             |
| `freeAll()` 置空指针          | ❌                    | ✅                             |
| `JOINT_NO` 定义在             | `.cpp`                | `.h`                           |
| `getRobotTorque()`            | ❌                    | ✅                             |

**`tests/observers.cpp` 是过时副本。** 它的 `freeAll()` 既不置空指针也不重置 `_nextIndex` —— **调用两次就是 double-free**。

**用 `matlab/observers.cpp`。**

---

## 9. 建议阅读顺序

按认知依赖排的，**不要按目录顺序读**。

### 第一遍 · 建立骨架（约 30 分钟）

1. **本文 §2 核心数学** —— 必须先推一遍 $\dot p = \tau - \beta + \tau_{ext}$，不然后面全是天书
2. [lib/external_observer.h](lib/external_observer.h) —— 全部 233 行，建立两条继承链的心智模型
3. [tests/double_link.h](tests/double_link.h) —— 看一个 `RobotDynamics` 的最小实现长什么样

### 第二遍 · 第一个观测器（最关键）

4. [lib/momentum_observer.h](lib/momentum_observer.h) —— 配合 §6.1 逐行对照。**这 98 行读透了，全库就通了一半**
5. [lib/momentum_observer_rnea.h](lib/momentum_observer_rnea.h) —— 和上面 diff 一下，只有两行不同。建立「RNEA 版 = 换动力学调用」的认知
6. [lib/robot_dynamics_rnea.cpp](lib/robot_dynamics_rnea.cpp) —— 配合 §5.2。**全库最精彩的 55 行**，务必看懂 `tranCqd` 怎么用 $\dot M = C + C^T$ 绕开转置
7. [tests/test.cpp](tests/test.cpp) —— 看观测器怎么被驱动，注意 `tau(0) -= 0.5` 的符号逻辑（§8.1）

### 第三遍 · 滤波器路线

8. [lib/iir_filter.h](lib/iir_filter.h) 的 `FilterF1` / `FilterF2` —— 配合 §5.3，**只看这两个**，其余三个是死代码
9. [lib/filtered_dyn_observer.h](lib/filtered_dyn_observer.h) —— 配合 §6.6。**理解 $\frac{\omega s}{s+\omega} = \omega - \frac{\omega^2}{s+\omega}$ 这个恒等式，就理解了 `FilterF2` 为什么长那样**
10. [lib/filtered_range_observer.h](lib/filtered_range_observer.h) —— 和上面 diff，注意「没有 tau」是**故意的**（§6.7）
11. [matlab/test.m:147-169](matlab/test.m#L147-L169) —— 看 `FRange` 的真实用途（自适应阈值）

### 第四遍 · 其余观测器（可按需跳读）

12. [lib/disturbance_observer.h](lib/disturbance_observer.h) —— §6.2，重点是 `L *= dt` 和隐式离散
13. [lib/sliding_mode_observer.h](lib/sliding_mode_observer.h) —— §6.3，重点是 super-twisting 和 `p` 的三次复用
14. [lib/kalman_filter.h](lib/kalman_filter.h) —— §5.4，注意两个 `step()` 重载的坑
15. [lib/disturbance_kalman_filter.h](lib/disturbance_kalman_filter.h) —— §6.4
16. [lib/kalman_filter_continous.cpp](lib/kalman_filter_continous.cpp) —— §5.5，Van Loan 方法（**注意 §10.1 的 bug**）
17. [lib/disturbance_kalman_filter_exp.h](lib/disturbance_kalman_filter_exp.h) —— §6.5

### 第五遍 · 集成（要做 MATLAB/ROS 对接才读）

18. [matlab/observers.cpp](matlab/observers.cpp) + [matlab/observers.h](matlab/observers.h) —— §8.2
19. [matlab/test.m](matlab/test.m) —— 完整的对比实验

### 可以跳过

- `lib/*_rnea.h` 里除 momentum 外的四个（模式完全一样）
- `iir_filter.h` 的后三个滤波器（死代码）
- `tests/observers.cpp`（过时副本）

---

## 10. 阅读时的坑（实测确认）

以下均已对照源码和 git 历史验证，不是猜测。

### 10.1 ⚠️ `KalmanFilterContinous` 读取未初始化内存

[lib/kalman_filter_continous.cpp:40-51](lib/kalman_filter_continous.cpp#L40-L51)：

```cpp
  // [A B;0 0]
  AB.resize(na+nb, na+nb);
  AB.setZero();                          // ← AB 清零了
  AB.block(0, 0, na, na) = a;
  AB.block(0, na, na, nb) = b;
  ABd.resize(na+nb, na+nb);
  // [-A Q; 0 A']
  AQ.resize(na+na, na+na);
  //                                     ← ★ AQ 没有 setZero()！
  AQ.block(0, 0, na, na) = -a;
  AQ.block(0, na, na, na) = Qd;
  AQ.block(na, na, na, na) = a.transpose();
  AQd.resize(na+na, na+na);
```

注释明明写着 `[-A Q; 0 A']`，**左下角那个 `0` 块从未被赋值**：

$$
\texttt{AQ} = \begin{bmatrix} -A & Q \\ \color{red}{\textbf{垃圾}} & A^T \end{bmatrix}
$$

Eigen 的 `resize()` 在尺寸改变时会重新分配且**不做初始化**（`AQ` 从构造时的 $1\times1$ 变成 $2n_a\times 2n_a$，必然重分配），所以 `AQ.block(na, 0, na, na)` 里是**未定义的垃圾值**。

紧邻的 `AB` 有 `setZero()`，`AQ` 没有 —— 一眼就能看出是遗漏。

**后果**：Van Loan 方法（§5.5）要求这个块严格为零，才能保证 $e^{\mathcal M}$ 是**块上三角**、从而

$$
Q_d = X_{22}^T X_{12}
$$

成立。垃圾值会通过 `exponential()` 的矩阵乘法**污染所有块**，导致 $Q_d$ 错误 → `DKalmanObserverExp` 结果不可信，且**每次运行结果可能不同**（取决于堆上的残留数据）。

**修复**：

```cpp
AQ.resize(na+na, na+na);
AQ.setZero();                 // ← 加这一行
AQ.block(0, 0, na, na) = -a;
...
```

**这个 bug 一直存在** —— `git show 675b5da^:lib/kalman_filter_continous.cpp` 确认早于 2024 年的代码风格重构，不是重构引入的。

**旁证**：`DKalmanObserverExp` 是 [matlab/test.m:106-126](matlab/test.m#L106-L126) 里**唯一被整段注释掉**的观测器。

### 10.2 ⚠️ `test.cpp` / `test_rnea.cpp` 的死分支编译不过

2024-12 的 `675b5da Update code style` 把 `Vector`/`Matrix` 改名成 `VectorJ`/`MatrixJ`。这次改名**漏掉了 `#ifdef` 里被预处理器排除的代码** —— 编译器根本看不到它们，所以没报错。


| 文件                                                     | 行                              | 内容                      | 所属分支 |
| -------------------------------------------------------- | ------------------------------- | ------------------------- | -------- |
| [tests/test.cpp:30](tests/test.cpp#L30)                  | `Vector k(2);`                  | `MomentumObserver`        | ❌       |
| [tests/test.cpp:47-50](tests/test.cpp#L47-L50)           | `Matrix S = Matrix::Zero(2,2);` | `DKalmanObserver`         | ❌       |
| [tests/test.cpp:56-59](tests/test.cpp#L56-L59)           | `Matrix S = ...`                | `DKalmanObserverExp`      | ❌       |
| [tests/test_rnea.cpp:39](tests/test_rnea.cpp#L39)        | `Vector T1(2), ...`             | `SlidingModeObserverRnea` | ❌       |
| [tests/test_rnea.cpp:47-50](tests/test_rnea.cpp#L47-L50) | `Matrix S = ...`                | `DKalmanObserverRnea`     | ❌       |

**当前恰好能编译**，因为启用的分支正好是干净的：

- `test.cpp` 开的是 `SlidingModeObserver`（用 `VectorJ` ✓）
- `test_rnea.cpp` 开的是 `MomentumObserverRnea`（用 `VectorJ` ✓）

**你一改注释切换观测器，就会撞上编译错误。** 把 `Vector` → `VectorJ`、`Matrix` → `MatrixJ` 即可。

同样过时的还有 [lib/README.md:5](lib/README.md#L5)：

```cpp
Vector getExternalTorque(Vector& q, Vector& dq, Vector& tau, double dt)
```

### 10.3 ⚠️ `cutOffHz` / `sampHz` 参数名是错的

```cpp
FDynObserver(RobotDynamics *rd, double cutOffHz, double sampHz);
FRangeObserver(RobotDynamics *rd, double cutOffHz, double sampHz, double k);
```

**两个名字都在撒谎**。证据链：

1. **底层实现的参数名是对的**：`FilterF1::update(double cutOff, double sampTime)` —— `sampTime` 不是 `sampHz`
2. **数学上**：`omega = tan(cutOff * sampTime * 0.5)` 要的是 $\tan\frac{\omega T}{2}$ → `cutOff` 必须是 **rad/s**，`sampTime` 必须是**秒**
3. **[matlab/test.m:131-132](matlab/test.m#L131-L132) 注释得明明白白**：
   ```matlab
   cutOff = 8;       % rad/s
   timeStep = 0.01;  % s
   ```
4. **[tests/test.cpp:65](tests/test.cpp#L65)**：`FDynObserver fd_observer(&robot, 8, TSTEP);`，而 `#define TSTEP 0.01`

**所以**：


| 参数名     | 声称 | 实际           | 差别      |
| ---------- | ---- | -------------- | --------- |
| `cutOffHz` | Hz   | **rad/s**      | $2\pi$ 倍 |
| `sampHz`   | Hz   | **秒（周期）** | 倒数关系  |

按名字理解传参会得到完全错误的滤波行为，而且**不会有任何报错**。

### 10.4 ⚠️ `KalmanFilter` 两个 `step()` 对 A 的假设相反

见 §5.4。`step(u,y)` 要**离散** A，`step(u,y,dt)` 要**连续** A。签名只差一个 `dt`，没有任何文档或断言提示。传错了静默出错。

附带：`step(u,y,dt)` 里 `A`、`B` 都按 `dt` 缩放了，但 **`Q` 没有** —— 意味着 `Q` 的含义是"每步"而非"每秒"，**控制周期一改必须重新整定 `Q`**。

### 10.5 变量复用严重损害可读性

这是读这个库**最大的实际障碍**。为了避免堆分配，同一个变量在一次函数调用里可能有三种含义：

[sliding_mode_observer.h](lib/sliding_mode_observer.h) 里的 `p`：

```cpp
p = dyn->getM(q) * qd;    // ① p = 动量
...
p = p_hat - p;            // ② p = 动量估计误差 p̃      ← "reuse p"
...
p = sigma;                // ③ p = 返回值（外部力矩）   ← "reuse to save result"
```

[momentum_observer.h](lib/momentum_observer.h) 里的 `torque`：

```cpp
torque = tau - dyn->getFriction(qd);   // ① τ - 摩擦
torque += r - beta;                    // ② 被积函数
torque(i) = r(i);                      // ③ 返回值
```

[robot_dynamics_rnea.cpp](lib/robot_dynamics_rnea.cpp) 里的 `_qext` / `_p0`：

```cpp
_qext.setZero(); _qext(i) = 1;    // _qext 当作单位向量 e_i（不是"扩展的 q"）
_p0 = rnea(q, _zero, _qext);      // _p0 当作矩阵的一列（不是"初始动量"）
```

**读的时候务必盯紧 `// reuse` 注释**，把每个变量的"当前含义"记在纸上。这些注释是作者留给你的路标，不是废话。

另外 `sum` 在 `MomentumObserver` 里身兼「积分累加器」和「$p(0)$」两职（§6.1），名字完全没体现。`KalmanFilter` 里的 `Y` 是新息协方差 $S$ 而参数 `y` 是观测量（§5.4）。

### 10.6 三个滤波器是死代码

`FilterButterworth`、`FilterLowPass`、`FilterHighPass` **没有任何观测器使用**（已 grep 全库确认）。而且 `FilterLowPass` 和 `FilterF1` 的系数、`filt()` 实现**完全相同**，纯属重复。

读 [lib/iir_filter.h](lib/iir_filter.h) 时只看 `FilterF1` 和 `FilterF2` 即可，后 100 行可以跳过。

### 10.7 `tests/observers.cpp` 是过时副本且有 double-free

见 §8.2 的对比表。`freeAll()` 既不置空指针也不重置 `_nextIndex`，**调用两次就 double-free**。以 [matlab/observers.cpp](matlab/observers.cpp) 为准。

### 10.8 `freeAll()` 的 switch 是多余的

见 §8.2。基类已有虚析构（[external_observer.h:150](lib/external_observer.h#L150)），`delete observer[i];` 一行就够。这个 switch 是 `7e4cf07 Add virtual destructors` 之前的遗留，加了虚析构后忘了简化。

**每加一个新观测器都必须记得改它，否则静默泄漏**（没有 `default` 分支）。

### 10.9 Kalman 观测器没有禁用拷贝

`DKalmanObserver` / `DKalmanObserverRnea` / `DKalmanObserverExp` 持有裸指针 `KalmanFilter *filter` 并在析构里 `delete`，但**没有禁用拷贝构造和拷贝赋值**。拷贝一个对象会导致 **double-free**。

改成 `std::unique_ptr<KalmanFilter>` 就自动安全（且自动禁用拷贝）。实践中不会这么用（都是 `new` 出来存指针），但这是个潜在陷阱。

### 10.10 其他

- **`GRAVITY = 9.81` 硬编码**在 [external_observer.h:17](lib/external_observer.h#L17)。RNEA 接口的 `g` 参数**默认是 0**，所以 `rnea(q, qd, q2d)` 不带重力 —— 读 `*_rnea.h` 时特别注意什么时候传了 `GRAVITY` 什么时候没传。
- **`tranCqd` 的注释 `// M'*qd - C*qd` 里的 `M'` 是 $\dot M = dM/dt$，不是转置**。在一个到处是 `.transpose()` 的文件里这么写，极易误读。
- **`DisturbanceObserver` 的三个参数最终塌缩成一个标量** $k = \frac{1}{2}(\xi + 2\beta\sigma)$，而它们的 Doxygen 注释全是 `@param sigma ...`（作者没写）。调参时建议固定 `sigma=1, xeta=0`，只调 `beta`。
- **`exponential()` 只有 5 项泰勒**（`#define EXP_TERMS 5`），不是精确矩阵指数。$\Delta t$ 大或 $A$ 特征值大时会不准，且无警告。
- **`DoubleLink::getFriction()` 返回零向量** —— 测试完全忽略摩擦，结果比真机乐观。
- **`DoubleLink::getC()` 的 `C(1,1)` 从未赋值**，靠构造函数的 `setZero()` 保持为 0。结果是对的，但很脆弱。
- **`KalmanFilter` 用简化协方差更新** $P = (I-KC)P$ 而非 Joseph 形式，长时间运行可能失去对称正定性。
- **`matlab/observers.cpp` 的封装层每次调用新建 4 个 `VectorJ`**，把库里精心设计的零分配努力全浪费了。照抄这个模式做实时节点要小心。
- **注释拼写错误**：`recursice`、`continous`（**贯穿文件名**！正确是 continuous）、`acceleraiton`、`estiamtion`、`@bried`、`Defauls`、`cinfiguration`、`Uncemment`、`expreimental`、`Reag`。搜代码时别按注释拼写搜。

---

## 11. 一句话总结全库

$$
\dot p = \tau - \beta + \tau_{ext}, \qquad \beta = G - C^T\dot q, \qquad p = M\dot q
$$

$$
\tau_{ext} = \dot p - \tau + \beta \quad\text{← 但 } \dot p \text{ 不能直接微分（噪声爆炸）}
$$

```
        ┌─────────────────────┼─────────────────────┐
        │                     │                     │
   【积分反馈】           【状态估计】            【滤波代数】
        │                     │                     │
  MomentumObserver      DKalmanObserver        FDynObserver
  ṙ = K(τ_ext - r)      X = [p; ω]             ωs/(s+ω) = ω - ω²/(s+ω)
  一阶低通                τ_ext = H·ω            → 微分被吸收
        │                     │                     │
  DisturbanceObserver   DKalmanObserverExp     FRangeObserver
  L = Y·M⁻¹ 加权          + Van Loan + R(q)      ×(±k) → 阈值包络带
  隐式离散，无条件稳定                                  │
        │                                             │
  SlidingModeObserver                          自适应阈值碰撞检测
  super-twisting                                      ↓
  有限时间收敛                              τ̂_ext 冲出包络 → 报警
```

```
  ═══════════════ 两种动力学接口，每个观测器两个版本 ═══════════════

  RobotDynamics                        RobotDynamicsRnea
  (显式给 M, C, G)                      (只有 rnea(q,q̇,q̈) 黑盒)
        │                                       │
   直接调用，快                          rnea(q,0,q̇)     → M·q̇   ← 速度塞进加速度位
                                        rnea(q,0,0,g)   → G
                                        rnea(q,q̇,0)     → C·q̇
                                        rnea(q,0,eᵢ)    → M 的第 i 列
                                        tranCqd(q,q̇)    → Cᵀq̇   ← N+2 次调用
                                              ↑
                                        用 Ṁ = C + Cᵀ 恒等式
                                        + 前向差分绕开转置
                                              │
                                        慢约 10 倍，但兼容
                                        Orocos / Pinocchio / KDL


  ═══════════════════ 与 dynamic_calibration 的接合 ═══════════════════

  dynamic_calibration                      ext_observer
  ────────────────────                     ─────────────
  estimate_drive_gains    ──→ K ──→  tau = K .* current    (test.m:37)
  estimate_dynamic_params ──→ π ──→  getM / getC / getG    (你的 RobotDynamics 实现)
                                              ↓
                                     getExternalTorque()
                                              ↓
                                           τ_ext
```
