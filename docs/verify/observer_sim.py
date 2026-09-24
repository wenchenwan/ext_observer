# -*- coding: utf-8 -*-
"""
ext_observer 算法的纯 Python 复现（无第三方依赖，Python 2.7 / 3.x 均可运行）。

用途：在没有 C++ 编译器和 Eigen 的环境里，逐行复现 lib/ 中观测器的离散算法，
在 tests/test.cpp 的同一仿真场景（2 连杆、1~2 s 注入 +0.5 N·m 外力）下验证：

  1. 各观测器能否恢复 +0.5 的外力方波；
  2. 动量观测器的离散稳定边界 K·dt < 2；
  3. 扰动观测器在测试参数下的离散极点（近似“无差拍”）；
  4. 两个卡尔曼观测器测试参数之间的等价关系。

运行：python observer_sim.py
说明：只复现 M/C/G 版本；摩擦为 0（与 tests/double_link.h 一致）。
"""
from __future__ import print_function
from math import sin, cos, tan, sqrt, tanh, fabs

# ---------------------------------------------------------------- 2 连杆模型
# 与 tests/double_link.h 完全相同的参数和公式
GRAV = 9.81
m1 = m2 = 1.0
l1 = 0.5
lc1 = lc2 = 0.25
I1, I2 = 0.3, 0.2


def M(q):
    a = m1*lc1*lc1 + m2*(l1*l1 + lc2*lc2 + 2*l1*lc2*cos(q[1])) + I1 + I2
    b = m2*(lc2*lc2 + l1*lc2*cos(q[1])) + I2
    c = m2*lc2*lc2 + I2
    return [[a, b], [b, c]]


def C(q, qd):
    h = -m2*l1*lc2*sin(q[1])
    return [[h*qd[1], h*(qd[0] + qd[1])], [-h*qd[0], 0.0]]


def G(q):
    c12 = cos(q[0] + q[1])
    return [(m1*lc1 + m2*l1)*GRAV*cos(q[0]) + m2*lc2*GRAV*c12, m2*lc2*GRAV*c12]


# ---------------------------------------------------------------- 2 维小工具
def mv(A, x):  return [A[0][0]*x[0] + A[0][1]*x[1], A[1][0]*x[0] + A[1][1]*x[1]]
def mtv(A, x): return [A[0][0]*x[0] + A[1][0]*x[1], A[0][1]*x[0] + A[1][1]*x[1]]  # A^T x
def add(a, b): return [a[0] + b[0], a[1] + b[1]]
def sub(a, b): return [a[0] - b[0], a[1] - b[1]]
def sc(s, a):  return [s*a[0], s*a[1]]


def inv(A):
    d = A[0][0]*A[1][1] - A[0][1]*A[1][0]
    return [[A[1][1]/d, -A[0][1]/d], [-A[1][0]/d, A[0][0]/d]]


# ---------------------------------------------------------------- 观测器
class MomentumObserver(object):          # lib/momentum_observer.h
    def __init__(self, k):
        self.k, self.run = k, False

    def step(self, q, qd, tau, dt):
        p = mv(M(q), qd)
        beta = sub(G(q), mtv(C(q, qd), qd))
        torque = tau[:]
        if self.run:
            torque = add(torque, sub(self.r, beta))
            self.sum = add(self.sum, sc(0.5*dt, add(torque, self.tprev)))
        else:
            torque = sub(torque, beta)
            self.r, self.sum, self.run = [0.0, 0.0], p[:], True
        self.tprev = torque
        p = sub(p, self.sum)
        self.r = [self.k[0]*p[0], self.k[1]*p[1]]
        return self.r[:]


class DisturbanceObserver(object):       # lib/disturbance_observer.h
    def __init__(self, sigma, xeta, beta):
        self.k, self.run = 0.5*(xeta + 2*beta*sigma), False

    def step(self, q, qd, tau, dt):
        Mi = inv(M(q))
        L = [[self.k*Mi[i][j]*dt for j in range(2)] for i in range(2)]
        p = sc(self.k, qd)
        if self.run:
            lft = [[1 + L[0][0], L[0][1]], [L[1][0], 1 + L[1][1]]]
            rht = add(self.z, mv(L, sub(sub(add(mv(C(q, qd), qd), G(q)), tau), p)))
            self.z = mv(inv(lft), rht)
        else:
            self.z, self.run = sc(-1, p), True
        return add(p, self.z)


class SlidingModeObserver(object):       # lib/sliding_mode_observer.h
    BIG = 50.0

    def __init__(self, T1, S1, T2, S2):
        self.T1, self.S1, self.T2, self.S2, self.run = T1, S1, T2, S2, False

    def step(self, q, qd, tau, dt):
        p = mv(M(q), qd)
        if not self.run:
            self.sig, self.ph, self.run = [0.0, 0.0], p[:], True
        pt = sub(self.ph, p)
        spp = [tanh(x*self.BIG) for x in pt]
        dph = add(sub(add(tau, mtv(C(q, qd), qd)), G(q)), self.sig)
        for i in range(2):
            dph[i] -= self.T2[i]*pt[i]
            dph[i] -= sqrt(fabs(pt[i]))*self.T1[i]*spp[i]
        ds = [-self.S1[i]*spp[i] - self.S2[i]*pt[i] for i in range(2)]
        out = self.sig[:]
        self.ph = add(self.ph, sc(dt, dph))
        self.sig = add(self.sig, sc(dt, ds))
        return out


class _F1(object):                       # lib/iir_filter.h: FilterF1
    def __init__(self, cut, T):
        self.cut = cut
        self.update(T)

    def update(self, T):
        w = tan(self.cut*T*0.5)
        self.k1, self.k2 = (1 - w)/(1 + w), w/(1 + w)

    def set(self, x):
        self.x1, self.y1 = x[:], x[:]

    def filt(self, x, dt):
        self.update(dt)
        self.y1 = add(sc(self.k1, self.y1), sc(self.k2, add(x, self.x1)))
        self.x1 = x[:]
        return self.y1[:]


class _F2(_F1):                          # lib/iir_filter.h: FilterF2
    def update(self, T):
        self.f2, self.w = 2.0/T, tan(self.cut*T*0.5)
        self.k1 = (1 - self.w)/(1 + self.w)
        self.k2 = -self.f2*self.w*self.w/(1 + self.w)

    def set(self, x):
        self.x1, self.y1 = x[:], sc(-self.f2*self.w, x)


class FDynObserver(object):              # lib/filtered_dyn_observer.h
    def __init__(self, cut, T):
        self.f1, self.f2, self.run = _F1(cut, T), _F2(cut, T), False

    def step(self, q, qd, tau, dt):
        p = mv(M(q), qd)
        rest = sub(sub(G(q), mtv(C(q, qd), qd)), tau)   # 摩擦为 0
        if self.run:
            res = add(self.f2.filt(p, dt), sc(self.f2.cut, p))
            return add(res, self.f1.filt(rest, dt))
        self.f2.set(p)
        self.f1.set(rest)
        self.run = True
        return [0.0, 0.0]


class DKalmanObserver(object):
    """lib/disturbance_kalman_filter.h + KalmanFilter::step(u, y, dt)，S=0, H=I。"""
    def __init__(self, Q, R):
        self.Q, self.R, self.run = Q, R, False

    def step(self, q, qd, tau, dt):
        n = 4
        p = mv(M(q), qd)
        u = add(sub(tau, G(q)), mtv(C(q, qd), qd))
        if not self.run:
            self.X = [p[0], p[1], 0.0, 0.0]
            self.P = [[0.0]*n for _ in range(n)]
            self.run = True
            return [0.0, 0.0]
        A = [[1, 0, dt, 0], [0, 1, 0, dt], [0, 0, 1, 0], [0, 0, 0, 1]]   # I + dt*[0 H; 0 S]
        X = [self.X[0] + dt*self.X[2] + dt*u[0], self.X[1] + dt*self.X[3] + dt*u[1],
             self.X[2], self.X[3]]
        AP = [[sum(A[i][k]*self.P[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
        P = [[sum(AP[i][k]*A[j][k] for k in range(n)) + self.Q[i][j] for j in range(n)]
             for i in range(n)]
        Yi = inv([[P[i][j] + self.R[i][j] for j in range(2)] for i in range(2)])
        K = [[sum(P[i][k]*Yi[k][j] for k in range(2)) for j in range(2)] for i in range(n)]
        nu = [p[0] - X[0], p[1] - X[1]]
        X = [X[i] + K[i][0]*nu[0] + K[i][1]*nu[1] for i in range(n)]
        P = [[P[i][j] - K[i][0]*P[0][j] - K[i][1]*P[1][j] for j in range(n)] for i in range(n)]
        self.X, self.P = X, P
        return [X[2], X[3]]


# ---------------------------------------------------------------- 仿真场景
def true_tau(q, qd, q2d, scale):
    """“真实”机器人的逆动力学：把观测器模型的 M、C、G 统一乘以 scale，用于制造模型误差。"""
    return sc(scale, add(add(mv(M(q), q2d), mv(C(q, qd), qd)), G(q)))


def simulate(obs, dt=0.01, t_end=3.0, ext=0.5, scale=1.0):
    """与 tests/test.cpp 相同：q = [sin 1.3t, sin 0.8t]，1<t<2 时 tau -= 0.5。
    scale != 1 时，真实机器人的动力学与观测器模型相差 (scale-1) 的比例。"""
    W1, W2 = 1.3, 0.8
    out, i, t = [], 0, 0.0
    while t < t_end:
        q = [sin(W1*t), sin(W2*t)]
        qd = [W1*cos(W1*t), W2*cos(W2*t)]
        q2d = [-W1*W1*sin(W1*t), -W2*W2*sin(W2*t)]
        tau = true_tau(q, qd, q2d, scale)
        if 1 < t < 2:
            tau = sub(tau, [ext, ext])
        out.append((t, obs.step(q, qd, tau, dt)))
        i += 1
        t = i*dt
    return out


def report(name, out):
    def at(tt):
        return min(out, key=lambda r: abs(r[0] - tt))[1]
    peak = max(max(abs(v) for v in r[1]) for r in out)
    print("%-26s t=0.5:[%7.3f %7.3f]  t=1.5:[%7.3f %7.3f]  t=2.5:[%7.3f %7.3f]  max|τ̂|=%.3g"
          % ((name,) + tuple(at(0.5)) + tuple(at(1.5)) + tuple(at(2.5)) + (peak,)))


if __name__ == "__main__":
    print("== 1. 外力恢复（期望：1~2 s 内约 +0.5，其余约 0）==")
    report("Momentum k=30", simulate(MomentumObserver([30, 30])))
    report("Disturbance 21/18/50", simulate(DisturbanceObserver(21, 18, 50)))
    S1, S2 = [15, 20], [10, 10]
    report("SlidingMode", simulate(SlidingModeObserver([2*sqrt(x) for x in S1], S1,
                                                       [2*sqrt(x) for x in S2], S2)))
    report("FDyn w=8", simulate(FDynObserver(8, 0.01)))
    Q = [[0.002, 0, 0, 0], [0, 0.002, 0, 0], [0, 0, 0.3, 0], [0, 0, 0, 0.3]]
    report("DKalman (Euler)", simulate(DKalmanObserver(Q, [[0.05, 0], [0, 0.05]])))

    print("\n== 2. 动量观测器稳定边界（理论：稳定 ⇔ K·dt < 2）==")
    for k in [50, 150, 190, 199, 201, 210, 250]:
        report("Momentum K·dt=%.2f" % (k*0.01), simulate(MomentumObserver([k, k])))

    print("\n== 3. 扰动观测器在测试参数下的离散极点 ==")
    k = 0.5*(18 + 2*50*21)
    Mq = M([0.0, 0.0])
    tr, det = Mq[0][0] + Mq[1][1], Mq[0][0]*Mq[1][1] - Mq[0][1]**2
    lam = [(tr + sqrt(tr*tr - 4*det))/2, (tr - sqrt(tr*tr - 4*det))/2]
    for lm in lam:
        rate = k/lm
        print("k=%.0f  eig(M)=%.4f  收敛速率 k/λ=%.0f 1/s  后向欧拉极点 1/(1+dt·k/λ)=%.4f"
              % (k, lm, rate, 1/(1 + 0.01*rate)))

    print("\n== 4. 两套卡尔曼参数的等价性（dt=0.01）==")
    dt = 0.01
    print("Exp 版 Q=0.2/30 → Q·dt = %.3f/%.3f（Euler 版为 0.002/0.3）" % (0.2*dt, 30*dt))
    print("Exp 版 R=5e-4  → R/dt = %.3f（Euler 版为 0.05，另乘 M·Mᵀ）" % (5e-4/dt))

    print("\n== 5. 动量观测器(K=w) 与 FDyn(w) 在连续时间上等价：对模型误差的反应也相同 ==")
    print("真实机器人动力学比模型大 10%，外力仍为 1~2 s 的 +0.5；从 t=1 s 起比较（避开两者初始化方式不同带来的启动瞬态）")
    w = 8.0
    for dt in [0.01, 0.001]:
        a = simulate(MomentumObserver([w, w]), dt=dt, scale=1.1)
        b = simulate(FDynObserver(w, dt), dt=dt, scale=1.1)

        def truth(t):
            return 0.5 if 1 < t < 2 else 0.0
        err_a = max(abs(r[1][j] - truth(r[0])) for r in a[int(1.0/dt):] for j in range(2))
        err_b = max(abs(r[1][j] - truth(r[0])) for r in b[int(1.0/dt):] for j in range(2))
        pairs = list(zip(a, b))[int(1.0/dt):]
        diff = max(abs(ra[1][j] - rb[1][j]) for ra, rb in pairs for j in range(2))
        print("dt=%-6g  MO 最大误差=%.4f  FDyn 最大误差=%.4f  两者之差最大=%.4f"
              % (dt, err_a, err_b, diff))
