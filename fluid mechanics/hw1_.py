import CFD 
import numpy as np
import numpy.linalg as npl
import matplotlib.pyplot as plt
import math


class CFD_hw1(CFD.CFD_base):
    def __init__(self, iterator_time= 1000):
        super(CFD_hw1, self).__init__(iterator_time)
    def initialize(self):
        self.gamma = 1.4 
        self.deltaX = 0.1
        start, end = 0, 3
        self.xs = np.linspace(start, end, (int)((end - start)//self.deltaX) + 2)
        self.sz = len(self.xs)
        xs = self.xs
        self.As = [1 + 2.2 * (x - 1.5) **2 for x in xs]
        self.ts = [0] # 时间
        self.rhos = [1 - 0.3146 * x for x in xs]
        self.Ts = [1 - 0.2314 * x for x in xs]
        self.Vs = [(0.1 + 1.09 * x) * np.sqrt(self.Ts[i]) for i, x in enumerate(xs)]
    def step(self):
        rhos = self.rhos
        Vs = self.Vs
        Ts = self.Ts
        As = self.As
        deltaX = self.deltaX
        gamma = self.gamma
        sz = self.sz
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
        self.ts += [self.ts[-1] + deltat] # ts = ts + [ts[-1] + deltat]
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
        self.rhos = rhoY + 0.5 * (rho_t + rho_t2) * deltat 
        self.Vs = VY + 0.5 * (V_t + V_t2) * deltat 
        self.Ts = TY + 0.5 * (T_t + T_t2) * deltat 
    def boundary_condition(self):
        self.rhos[0] = 1 
        self.Ts[0] = 1
        self.Vs[0] = 2 * self.Vs[1] - self.Vs[2] 
        self.rhos[-1] = 2 * self.rhos[-2] - self.rhos[-3]
        self.Ts[-1] = 2 * self.Ts[-2] - self.Ts[-3] 
        self.Vs[-1] = 2 * self.Vs[-2] - self.Vs[-3]
    def solve_return_value(self):
        rev = (self.rhos, self.Vs, self.Ts, self.ts, self.xs)
        return rev

if __name__ == '__main__':
    cfd_hw1 = CFD_hw1(1000)
    rhos, Vs, Ts, ts, xs = cfd_hw1.solve()
    print(rhos)
    print(Vs)
    plt.plot(xs, rhos)
    plt.savefig('1.png')

    cfd_hw12 = CFD.CFD_MacCormack()
    rhos, Vs, Ts, ts, xs = cfd_hw12.solve()
    print(rhos)
    print(Vs)
    plt.plot(xs, rhos)
    plt.savefig('2.png')
