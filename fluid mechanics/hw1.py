import numpy as np
import numpy.linalg as npl
import matplotlib.pyplot as plt

gamma = 1.4 
deltaX = 0.1
start, end = 0, 3
xs = np.linspace(start, end, (int)((end - start)//deltaX) + 2)
print(xs)

print(xs[0])
print(xs[-1])
print(xs[-2])
print(xs[-2])

rhos = [1 - 0.3146 * x for x in xs]
Ts = [1 - 0.2314 * x for x in xs] # 温度
Vs = [(0.1 + 1.09 * x) * np.sqrt(Ts[i]) for i, x in enumerate(xs)]
print(rhos)

print([1, 3, 41])
# [1, 1, 1] + [0, 1, 2] = [1, 2, 3]
print([1, 1, 1] + [0, 1, 2])
n1 = np.array([1, 1, 1])
n2 = np.array([0, 1, 2])
print(np.array([1, 1, 1]) + np.array([0, 1, 2]))

sz = len(xs)
As = [1 + 2.2 * (x - 1.5) **2 for x in xs]
print(As)

ts = [0] # 时间
rhos = [1 - 0.3146 * x for x in xs]
Ts = [1 - 0.2314 * x for x in xs]
Vs = [(0.1 + 1.09 * x) * np.sqrt(Ts[i]) for i, x in enumerate(xs)]
rho_mid = [0.5]
V_mid = [Vs[16]]
T_mid = [Ts[16]]
rho_t_mid = [0]
rho_t_mid2 = [0]

for it in range(1000):
    rho_t = np.array([-rhos[i] * (Vs[i + 1] - Vs[i]) / deltaX - Vs[i] * (rhos[i + 1] - rhos[i]) / deltaX - rhos[i] * Vs[i] * (np.log(As[i + 1]) - np.log(As[i])) / deltaX for i in range(sz - 1)])
    V_t = np.array([- Vs[i] * (Vs[i + 1] - Vs[i]) / deltaX - 1 / gamma * ((Ts[i + 1] - Ts[i]) / deltaX + Ts[i] / rhos[i] * (rhos[i + 1] - rhos[i]) / deltaX) for i in range(sz - 1)])
    T_t = np.array([- Vs[i] * (Ts[i + 1] - Ts[i]) / deltaX - (gamma - 1) * Ts[i] * ((Vs[i + 1] - Vs[i]) / deltaX + Vs[i] * (np.log(As[i + 1]) - np.log(As[i])) / deltaX) for i in range(sz - 1)])
    deltats = [deltaX / (Vs[i] + np.sqrt(Ts[i])) for i in range(sz - 1)]
    deltat = 0.9 * np.min(deltats)
    # print(deltats)
    # print(Ts)
    # print(Vs)
    # print(rhos)
    # print("-------------------------------")
    ts += [ts[-1] + deltat] # ts = ts + [ts[-1] + deltat]
    rho_t2 = np.append(rho_t, rho_t[-1])
    V_t2 = np.append(V_t, V_t[-1])
    T_t2 = np.append(T_t, T_t[-1])

    rhoY = rhos 
    VY = Vs 
    TY = Ts 
    rhos = rhos + rho_t2 * deltat 
    Vs = Vs + V_t2 * deltat 
    Ts = Ts + T_t2 * deltat 

    rho_t = np.array([-rhos[i] * (Vs[i] - Vs[i - 1]) / deltaX - Vs[i] * (rhos[i] - rhos[i - 1]) / deltaX - rhos[i] * Vs[i] * (np.log(As[i]) - np.log(As[i - 1])) / deltaX for i in range(1, sz)])
    V_t = np.array([- Vs[i] * (Vs[i] - Vs[i - 1]) / deltaX - 1 / gamma * ((Ts[i] - Ts[i - 1]) / deltaX + Ts[i] / rhos[i] * (rhos[i] - rhos[i - 1]) / deltaX) for i in range(1, sz)])
    T_t = np.array([- Vs[i] * (Ts[i] - Ts[i - 1]) / deltaX - (gamma - 1) * Ts[i] * ((Vs[i] - Vs[i - 1]) / deltaX + Vs[i] * (np.log(As[i]) - np.log(As[i - 1])) / deltaX) for i in range(1, sz)])

    rho_t = np.hstack([rho_t[0], rho_t])# [1, 2, 0.7, ...] -> [1, 1, 2, 0.7, ...]
    V_t = np.hstack([V_t[0], V_t])
    T_t = np.hstack([T_t[0], T_t])

    # 计算时间偏导数的平均值，计算矫正值
    rhos = rhoY + 0.5 * (rho_t + rho_t2) * deltat 
    Vs = VY + 0.5 * (V_t + V_t2) * deltat 
    Ts = TY + 0.5 * (T_t + T_t2) * deltat 

    # 边界条件
    rhos[0] = 1 
    Ts[0] = 1
    Vs[0] = 2 * Vs[1] - Vs[2] 
    rhos[-1] = 2 * rhos[-2] - rhos[-3]
    Ts[-1] = 2 * Ts[-2] - Ts[-3] 
    Vs[-1] = 2 * Vs[-2] - Vs[-3]

    # 记录喉道数据
    rho_mid += [rhos[16]]
    V_mid += [Vs[16]]
    T_mid += [Ts[16]]
    rho_t_mid += [rho_t[16]]
    rho_t_mid2 += [rho_t2[16]]

print(ts)
print(rhos)

 

plt.plot(ts, rho_mid)
plt.plot(ts, V_mid)
plt.plot(ts, T_mid)
plt.legend(['rho', 'V', 'T'])

plt.plot(ts, rho_t_mid)
plt.plot(ts, rho_t_mid2)
plt.plot(ts, np.array(rho_t_mid2) - np.array(rho_t_mid))

plt.plot(ts[100:110], rho_t_mid[100:110])
plt.plot(ts[100:110], rho_t_mid2[100:110])
plt.plot(ts[100:110], (np.array(rho_t_mid2) - np.array(rho_t_mid))[100:110])

plt.plot(xs, rhos)
cs = [np.sqrt(T) for T in Ts]
rhos_ = [(1 + 0.2 * (Vs[i] / cs[i]) ** 2) ** ( -1 / 0.4) for i in range(len(rhos))]
plt.plot(xs, rhos_)

plt.plot(xs, Vs)

plt.plot(xs, Ts)