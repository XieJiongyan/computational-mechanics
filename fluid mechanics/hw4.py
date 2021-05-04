import CFD 
import numpy as np
import numpy.linalg as npl
import matplotlib.pyplot as plt
import math

class CFD_hw4(CFD.CFD_base):
    def __init__(self, iterator_time=1000):
        super().__init__(iterator_time=iterator_time) 
    def initialize(self):
        self.h = 0.25
        self.phi = np.zeros([15, 9]) #15个横向点，9个纵向点
        self.phi_ = self.phi #上一个时刻的phi
        self.lattice_type = np.zeros([15, 9]) # 0: inside, 1: not exactly on boundary, 2: just on boundary, 3: outside
        for i in range(14, -1, -1):
            for j in range(9):
                if (3.5 - i * self.h) ** 2 + (j * self.h) ** 2 < 1:
                    self.lattice_type[i, j], self.lattice_type[i - 1, j], self.lattice_type[i, j + 1] = 3, 1, 1
                elif (3.5 - i * self.h) ** 2 + (j * self.h) ** 2 == 1:
                    self.lattice_type[i, j] = 2 
    
    def boundary_condition(self):
        for i in range(0, 15):
            if self.lattice_type[i, 0] == 0:
                self.phi[i, 0] = 0
            if self.lattice_type[i, -1] == 0:
                self.phi[i, -1] = 2 
        for j in range(0, 8):
            if (self.lattice_type[-1, j]) == 0:
                self.phi[-1, j] = 0.25 * (2 * self.phi_[-2, j] +  self.phi_[-1, j - 1] + self.phi_[-1, j + 1])
            self.phi[0, j] = j * self.h 
        h = self.h 
        for i in range(15):
            for j in range(9):
                if self.lattice_type[i, j] == 2:
                    self.phi[i, j] = 0
                elif self.lattice_type[i, j] == 1: 
                    a = 15 * self.h - math.sqrt(1 - (j * self.h) ** 2) 
                    b = j * h - math.sqrt(1 - ((14 - i) * self.h) ** 2) 
                    psi1, psi2 = 0, 0
                    if a > self.h:
                        psi1 = self.phi_[i + 1, j] 
                    if b > self.h:
                        psi2 = self.phi_[i, j - 1] 
                    self.phi[i, j] = (self.phi_[i - 1, j]  / h / (a + h) + self.phi_[i, j + 1] / h / (b + h) + psi1 / a / (a + h) + psi2 / b / (b + h)) / ( 1 / a / h + 1 / b / h)
                    
        self.phi_ = self.phi 
    def step(self):
        for i in range(1, 14):
            for j in range(1, 8):
                if self.lattice_type[i, j] == 0:
                    self.phi[i, j] = 0.25 * (self.phi_[i - 1, j] + self.phi_[i, j - 1] + self.phi_[i + 1, j] + self.phi_[i, j + 1] )
    def solve_return_value(self):
        return self.phi
    
                    
        
if __name__ == '__main__':
    '''a'''
    cfd = CFD_hw4() 
    cfd.solve()
    print(cfd.lattice_type) 
    print(np.around(cfd.phi, 2))
    plt.contour(cfd.phi.T, np.linspace(0.5, 2, 7)) 
    plt.savefig('fig/CFD/hw4.png')

    