import CFD 
import numpy as np
import numpy.linalg as npl
import matplotlib.pyplot as plt
import math

class CFD_hw4(CFD.CFD_base):
    def __init__(self, iterator_time=1000):
        super().__init__(iterator_time=iterator_time) 
    def initialize(self):
        length, width = 3.5, 2
        self.h = 0.10
        self.ls, self.ws = int(length // self.h + 1), int(width // self.h + 1)
        self.phi = np.zeros([self.ls, self.ws]) #15个横向点，9个纵向点
        self.phi_ = self.phi #上一个时刻的phi
        self.lattice_type = np.zeros([self.ls, self.ws]) # 0: inside, 1: not exactly on boundary, 2: just on boundary, 3: outside
        for i in range(self.ls - 1, -1, -1):
            for j in range(self.ws):
                if (3.5 - i * self.h) ** 2 + (j * self.h) ** 2 < 1 - 1e-5:
                    self.lattice_type[i, j], self.lattice_type[i - 1, j], self.lattice_type[i, j + 1] = 3, 1, 1
                elif abs((3.5 - i * self.h) ** 2 + (j * self.h) ** 2  - 1) < 1e-5:
                    self.lattice_type[i, j] = 2 
    
    def boundary_condition(self):
        for i in range(0, self.ls):
            if self.lattice_type[i, 0] == 0:
                self.phi[i, 0] = 0
            if self.lattice_type[i, -1] == 0:
                self.phi[i, -1] = 2 
        for j in range(0, self.ws - 1):
            if (self.lattice_type[-1, j]) == 0:
                self.phi[-1, j] = 0.25 * (2 * self.phi_[-2, j] +  self.phi_[-1, j - 1] + self.phi_[-1, j + 1])
            self.phi[0, j] = j * self.h 
        h = self.h 
        for i in range(self.ls):
            for j in range(self.ws):
                if self.lattice_type[i, j] == 2:
                    self.phi[i, j] = 0
                elif self.lattice_type[i, j] == 1: 
                    if i == self.ls - 1:
                        continue 
                    a = self.ls * self.h - math.sqrt(1 - (j * self.h) ** 2) 
                    b = j * h - math.sqrt(1 - ((self.ls - 1 - i) * self.h) ** 2) 
                    psi1, psi2 = 0, 0
                    if a > self.h:
                        psi1 = self.phi_[i + 1, j] 
                    if b > self.h:
                        psi2 = self.phi_[i, j - 1] 
                    self.phi[i, j] = (self.phi_[i - 1, j]  / h / (a + h) + self.phi_[i, j + 1] / h / (b + h) + psi1 / a / (a + h) + psi2 / b / (b + h)) / ( 1 / a / h + 1 / b / h)
        self.phi_ = self.phi 
        print("step")
        print(np.around(cfd.phi, 2))
        print(np.around(cfd.phi_, 2))
    def step(self):
        for i in range(1, self.ls - 1):
            for j in range(1, self.ws - 1):
                if self.lattice_type[i, j] == 0:
                    self.phi[i, j] = 0.25 * (self.phi_[i - 1, j] + self.phi_[i, j - 1] + self.phi_[i + 1, j] + self.phi_[i, j + 1] )
    def solve_return_value(self):
        return self.phi
    
                    
        
if __name__ == '__main__':
    '''a'''
    cfd = CFD_hw4(2) 
    cfd.solve()
    print(cfd.lattice_type) 
    # print(np.around(cfd.phi, 2))
    plt.contour(cfd.phi.T, np.linspace(0, 2, 8)) 
    plt.savefig('fig/CFD/hw4.png')

    