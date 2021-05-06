import CFD 
import numpy as np
import numpy.linalg as npl
import matplotlib.pyplot as plt
import math

class CFD_hw5(CFD.CFD_base):
    def __init__(self, iterator_time=1000):
        super().__init__(iterator_time=iterator_time) 
    def initialize(self):
        self.m = 11 
        self.h = 1 / self.m 
        self.psi = np.zeros([self.m, self.m]) 
        self.zeta = np.zeros([self.m, self.m]) 

    
    def boundary_condition(self):
        for i in range(0, self.m):
            self.psi[i, 0], self.psi[i, -1], self.psi[0, i], self.psi[-1, i] = 0, 0, 0, 0 
        for i in range(1, self.m - 1): 
            self.zeta[i, 0] = 1 / self.h ** 2 * (-self.psi[i - 1, 1] + 8 / 3 * self.psi[i, 1] - self.psi[i + 1, 1] - 2 / 3 * self.psi[i, 3]) + 2 / 3 / self.h
            self.zeta[i, -1] = 0 
            self.zeta[0, i] = 1 / self.h ** 2 * (-self.psi[1, i - 1] + 8 / 3 * self.psi[1, i] - self.psi[1, i + 1] - 2 / 3 * self.psi[2, i])
            self.zeta[-1, i] = 1 / self.h ** 2 * (-self.psi[-2, i - 1] + 8 / 3 * self.psi[-2, i] - self.psi[-2, i + 1] - 2 / 3 * self.psi[-3, i])
        self.zeta[0, 0], self.zeta[0, -1], self.zeta[-1, 0], self.zeta[-1, -1] = 0.5 * (self.zeta[0, 1] + self.zeta[1, 0]), 0.5 * (self.zeta[0, -2] + self.zeta[1, -1]), 0.5 * (self.zeta[-2, 0] + self.zeta[-1, 1]), 0.5 * (self.zeta[-2, -1] + self.zeta[-1, -2])
                    
    def step(self):
        self.get_zeta() 
        self.get_psi()
    def get_zeta(self): 
        for t in range(1000):
            zeta = self.zeta 
            for i in range(1, self.m - 1):
                for j in range(1, self.m - 1):
                    self.zeta[i, j] = 0.25 * (zeta[i - 1, j] + zeta[i + 1, j] + zeta[i, j - 1] + zeta[i, j + 1]) 
            if npl.norm(zeta - self.zeta) < 1e-3:
                break 
        else: 
            print("over 1000 iteration") 

    def get_psi(self):
        for t in range(1000):
            psi = self.psi 
            for i in range(1, self.m - 1): 
                for j in range(1, self.m - 1): 
                    self.psi[i, j] = 0.25 * (self.zeta[i, j] * self.h ** 2 + self.psi[i - 1, j] + self.psi[i + 1, j] + self.psi[i, j - 1] + self.psi[i, j + 1])
            if npl.norm(psi - self.psi) < 1e-3: 
                break 
        else:
            print("over 1000 iterations")
    def solve_return_value(self):
        return self.psi, self.zeta

if __name__ == '__main__':
    hw5 = CFD_hw5()
    psi, zeta = hw5.solve() 
    print(np.round(psi, 2))
    plt.contour(psi.T, np.linspace(0, 0.1, 5)) 
    plt.savefig('fig/CFD/hw5.png')
