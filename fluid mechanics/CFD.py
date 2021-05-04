import numpy as np
import numpy.linalg as npl
import matplotlib.pyplot as plt
import math

class CFD_base():
    def __init__(self, iterator_time= 1000):
        '''init'''
        self.iterator_time = iterator_time
    def initialize(self):
        '''initialize'''
    def boundary_condition(self):
        '''boundary_condition'''
    def stop_condition(self):
        return False
    def step(self):
        '''step'''        
    
    # 返回值
    def solve_return_value():
        '''..'''
        return None
    def solve(self):
        self.initialize()
        for it in range(self.iterator_time):
            self.step()
            self.boundary_condition()
            if (self.stop_condition()):
                break
        return self.solve_return_value()

class CFD_MacCormack(CFD_base):
    def __init__(self, iterator_time= 1000, deltaX= 0.1):
        super(CFD_MacCormack, self).__init__(iterator_time)
        self.deltaX = deltaX
        self.Kelang = 0.9
    def initial_condition(self):
        xs = self.xs
        self.rhos = [1 - 0.3146 * x for x in xs]
        self.Ts = [1 - 0.2314 * x for x in xs]
        self.Vs = [(0.1 + 1.09 * x) * np.sqrt(self.Ts[i]) for i, x in enumerate(xs)]
    def setAs(self):
        self.As = [1 + 2.2 * (x - 1.5) **2 for x in self.xs]
    # @brief  画网格，设置gamma
    def initialize(self):
        self.gamma = 1.4 
        start, end = 0, 3
        self.xs = np.linspace(start, end, (int)((end - start)//self.deltaX) + 2)
        self.sz = len(self.xs)
        self.ts = [0] # 时间
        self.setAs()
        self.initial_condition()
    def calculate_step(self, i, x, y):
        rhos = self.rhos
        Vs = self.Vs
        Ts = self.Ts
        As = self.As
        deltaX = self.deltaX
        gamma = self.gamma
        sz = self.sz
        rho_t = -rhos[i] * (Vs[x] - Vs[y]) / deltaX - Vs[i] * (rhos[x] - rhos[y]) / deltaX - rhos[i] * Vs[i] * (np.log(As[x]) - np.log(As[y])) / deltaX
        V_t = - Vs[i] * (Vs[x] - Vs[y]) / deltaX - 1 / gamma * ((Ts[x] - Ts[y]) / deltaX + Ts[i] / rhos[i] * (rhos[x] - rhos[y]) / deltaX) 
        T_t = - Vs[i] * (Ts[x] - Ts[y]) / deltaX - (gamma - 1) * Ts[i] * ((Vs[x] - Vs[y]) / deltaX + Vs[i] * (np.log(As[x]) - np.log(As[y])) / deltaX) 
        return rho_t, V_t, T_t 
    def forward_difference(self):
        rho_t = np.zeros([self.sz])
        V_t = np.zeros([self.sz])
        T_t = np.zeros([self.sz])
        for i in range(self.sz - 1):
           rho_t[i], V_t[i], T_t[i] = self.calculate_step(i, i + 1, i)
        return rho_t, V_t, T_t 
    def backward_difference(self):
        rho_t = np.zeros([self.sz])
        V_t = np.zeros([self.sz])
        T_t = np.zeros([self.sz])
        for i in range(1, self.sz):
           rho_t[i], V_t[i], T_t[i] = self.calculate_step(i, i, i - 1)
        return rho_t, V_t, T_t 
    def step_t(self):
        deltats = [self.deltaX / (self.Vs[i] + np.sqrt(self.Ts[i])) for i in range(self.sz - 1)]
        deltat = self.Kelang * np.min(deltats)
        self.ts += [self.ts[-1] + deltat] # ts = ts + [ts[-1] + deltat]
        return deltat
    def step(self):
        deltat = self.step_t()
        rho_t, V_t, T_t = self.forward_difference()

        # 保存现场的rhos, Vs, Ts, 并更新rhos, Vs, Ts
        rhoY = self.rhos 
        VY = self.Vs 
        TY = self.Ts 
        self.rhos = self.rhos + rho_t * deltat 
        self.Vs = self.Vs + V_t * deltat 
        self.Ts = self.Ts + T_t * deltat 

        rho_t2, V_t2, T_t2 = self.backward_difference()

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

class CFD_Mac_conservation(CFD_MacCormack):
    def __init__(self, iterator_time=1000, deltaX=0.1):
        super().__init__(iterator_time=iterator_time, deltaX=deltaX)
    def initial_condition(self):
        super().initial_condition()
        self.U1 = [0] * self.sz 
        self.U2 = [0] * self.sz 
        self.U3 = [0] * self.sz 
        xs = self.xs
        self.rhos = [1 if x <= 0.5 else 1 - 0.366 * (x - 0.5) if x <= 0.5 else 0.634 - 0.3879 * (x - 1.5) for x in xs]
        self.Ts = [1 if x <= 0.5 else 1.0 - 0.167 * (x - 0.5) if x <= 1.5 else 0.833 - 0.3507 * (x - 1.5) for x in xs]
        self.Vs = [0.59 / self.rhos[i] / self.As[i] for i in range(self.sz)]
        self.U1 = [self.rhos[i] * self.As[i] for i in range(self.sz)]
        self.U2 = [self.rhos[i] * self.As[i] * self.Vs[i] for i in range(self.sz)]
        self.U3 = [self.rhos[i] * (self.Ts[i] / (self.gamma - 1) + self.gamma / 2 * self.Vs[i] ** 2) * self.As[i] for i in range(self.sz)]

    def calculate_step(self, i, x, y):
        U1_t = - (self.F1[x] - self.F1[y]) / self.deltaX
        U2_t = - (self.F2[x] - self.F2[y]) / self.deltaX
        U3_t = - (self.F3[x] - self.F3[y]) / self.deltaX
        return U1_t, U2_t, U3_t
    # @usage self.getF_Ut()
    def getF_Ut(self):
        U1, U2, U3 = self.U1, self.U2, self.U3 
        self.F1 = U2 
        self.F2 = [U2[i] ** 2 / U1[i] + (self.gamma - 1) / self.gamma * (U3[i] - self.gamma / 2 * (U2[i] / U1[i]) ** 2) for i in range(self.sz)]
        self.F3 = [self.gamma * U2[i] * U3[i] / U1[i] - self.gamma * (self.gamma - 1) / 2 * U2[i] ** 3 / U1[i] ** 2 for i in range(self.sz)]
        self.J2 = [1 / self.gamma * self.rhos[i] * self.Ts[i] * (self.As[i + 1] - self.As[i]) / self.deltaX for i in range(self.sz - 1)]
    def update_origin_variables(self):
        U1, U2, U3 = self.U1, self.U2, self.U3 
        self.rhos = [U1[i] / self.As[i] for i in range(self.sz)]
        self.Vs = [U1[i] / U2[i] for i in range(self.sz)]
        self.Ts = [(self.gamma - 1) * (U3[i] / U1[i] - self.gamma / 2 * (U2[i] / U1[i]) ** 2) for i in range(self.sz)]
        self.ps = [self.rhos[i] * self.Ts[i] for i in range(self.sz)]
    def step(self):
        deltat = self.step_t()
        self.getF_Ut()
        U1_t, U2_t, U3_t = self.forward_difference()
        self.update_origin_variables()

        # 保存现场的rhos, Vs, Ts, 并更新rhos, Vs, Ts
        U1 = self.U1 
        U2 = self.U2 
        U3 = self.U3 
        self.U1 = [self.U1[i] + U1_t[i] * deltat for i in range(self.sz)] 
        self.U2 = [self.U2[i] + U2_t[i] * deltat for i in range(self.sz)] 
        self.U3 = [self.U3[i] + U3_t[i] * deltat for i in range(self.sz)] 

        self.getF_Ut()
        U1_t2, U2_t2, U3_t2 = self.forward_difference()

        # 计算时间偏导数的平均值，计算矫正值
        self.U1 = [U1[i] + 0.5 * (U1_t[i] + U1_t2[i]) * deltat for i in range(self.sz)]
        self.U2 = [U2[i] + 0.5 * (U2_t[i] + U2_t2[i]) * deltat for i in range(self.sz)]
        self.U3 = [U3[i] + 0.5 * (U3_t[i] + U3_t2[i]) * deltat for i in range(self.sz)]
        self.update_origin_variables()
    def boundary_condition(self):
        self.U1[0] = self.As[0]
        self.U2[0] = 2 * self.U2[1] - self.U2[2]
        self.U3[0] = self.U1[0] * (self.Ts[0] / (self.gamma - 1) + self.gamma / 2 * self.Vs[0] ** 2 )

        self.U1[-1] = 2 * self.U1[-2] - self.U1[-3]
        self.U2[-1] = 2 * self.U2[-2] - self.U2[-3]
        self.U2[-1] = 2 * self.U1[-2] - self.U2[-3]
