import numpy as np
import numpy.linalg as npl
import matplotlib.pyplot as plt

gamma = 1.4 
deltaX = 0.1
start, end = 0, 3
xs = np.linspace(start, end, (int)((end - start)//deltaX) + 2)
print(xs)
sz = xs.shape[0]
As = np.array([0.0] * sz)
mid1_5 = int(sz // 2 + 1) # x = 1.5时对应的节点
print(mid1_5)
for i in range(mid1_5):
    As[i] = 1 + 2.2 * (xs[i] - 1.5) ** 2
for i in range(mid1_5, sz):
    As[i] = 1 + 0.2223 * (xs[i] - 1.5) ** 2
print(As)
pN = (1 - 0.023 * 3) * (1 - 0.000933 * 3)
pe = 0.93

ts = [0] # 时间
rhos = [1 - 0.023 * x for x in xs]
print(1 - 0.023 * 0.1)
Ts = [1 - 0.009333 * x for x in xs]
Vs = [(0.05 + 0.11 * x) * np.sqrt(1 - 0.2314 * x) for i, x in enumerate(xs)]

print(0.05 + 0.11 * 0.1 * np.sqrt(1 - 0.2314 * 0.1))
rho_mid = [rhos[16]]
V_mid = [Vs[16]]
T_mid = [Ts[16]]
rho_t_mid = [0]
rho_t_mid2 = [0]

for it in range(5000):
    print('rhos: ', rhos)
    print('Ts: ', Ts)
    print('Vs: ', Vs)
    rho_t = np.array([-rhos[i] * (Vs[i + 1] - Vs[i]) / deltaX - Vs[i] * (rhos[i + 1] - rhos[i]) / deltaX - rhos[i] * Vs[i] * (np.log(As[i + 1]) - np.log(As[i])) / deltaX for i in range(sz - 1)])
    V_t = np.array([- Vs[i] * (Vs[i + 1] - Vs[i]) / deltaX - 1 / gamma * ((Ts[i + 1] - Ts[i]) / deltaX + Ts[i] / rhos[i] * (rhos[i + 1] - rhos[i]) / deltaX) for i in range(sz - 1)])
    T_t = np.array([- Vs[i] * (Ts[i + 1] - Ts[i]) / deltaX - (gamma - 1) * Ts[i] * ((Vs[i + 1] - Vs[i]) / deltaX + Vs[i] * (np.log(As[i + 1]) - np.log(As[i])) / deltaX) for i in range(sz - 1)])
    print(rho_t)
    print(V_t)
    print(T_t)
    print('---------')
    deltats = [deltaX / (Vs[i] + np.sqrt(Ts[i])) for i in range(1, sz - 1)]
    deltat = 0.5 * np.min(deltats)
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
    Vs[-1] = 2 * Vs[-2] - Vs[-3]
    # Ts[-1] = 2 * Ts[-2] - Ts[-3] 
    # rhos[-1] = 0.93 / Ts[-1]
    rhos[-1] = 2 * rhos[-2] - rhos[-3]
    Ts[-1] = pe / rhos[-1]
    print(rhos[-1], Ts[-1], rhos[-1] * Ts[-1])

    print(rhos)
    print(Vs)
    print(Ts)
    # 记录喉道数据
    rho_mid += [rhos[16]]
    V_mid += [Vs[16]]
    T_mid += [Ts[16]]
    rho_t_mid += [rho_t[16]]
    rho_t_mid2 += [rho_t2[16]]

import math
Mas = [Vs[i] / math.sqrt(Ts[i]) for i in range(sz)]
print(Mas)
Mae=np.sqrt((pe**((-gamma+1)/gamma)-1)*2/(gamma-1))
A_= As[-1] / np.sqrt(Mae**(-2)*(2/(gamma+1)*(1+(gamma-1)/2*Mae**2))**((gamma+1)/(gamma-1)))
print('A_ == ', A_)
from scipy.optimize import fsolve 
Ma_exact = []
for i in range(sz):
    x=deltaX*i
    thisA=As[i]
    def func(Ma:int):
        return (thisA/A_)**2-Ma**(-2)*(2/(gamma+1)*(1+(gamma-1)/2*Ma**2))**((gamma+1)/(gamma-1))
    # print(fsolve(func,Mas[i]))
    Ma_exact.append(fsolve(func,Mas[i])[0])
print(Ma_exact)

plt.plot(ts, rho_mid)
plt.plot(ts, V_mid)
plt.plot(ts, T_mid)
plt.legend(['rho', 'V', 'T'])

plt.plot(xs, rhos)
plt.plot(xs, Vs, )
plt.plot(xs, Ts)
plt.plot(xs, Mas, 'o')
plt.plot(xs, Ma_exact, 'o')
plt.plot(xs, [Ts[i] * rhos[i] for i in range(sz)])
plt.legend([r'$\rho$', r'$v$', r'$T$', r'Calculate $Ma$', r'Exact $Ma$', r'$c$'])
plt.show()