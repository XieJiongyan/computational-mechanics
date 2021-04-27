import CFD 
import numpy as np
import numpy.linalg as npl
import matplotlib.pyplot as plt
import math

class CFD_hw2(CFD.CFD_MacCormack):
    def __init__(self, iterator_time= 1000):
        super(CFD_hw2, self).__init__(iterator_time)
        self.Kelang = 0.5
        self.pe = 0.93
    def setAs(self):
        sz = self.sz 
        self.As = np.array([0.0] * sz)
        mid1_5 = int(sz // 2 + 1) # x = 1.5时对应的节点
        print(mid1_5)
        for i in range(mid1_5):
            self.As[i] = 1 + 2.2 * (self.xs[i] - 1.5) ** 2
        for i in range(mid1_5, sz):
            self.As[i] = 1 + 0.2223 * (self.xs[i] - 1.5) ** 2
    def initial_condition(self):
        self.rhos = [1 - 0.023 * x for x in self.xs]
        self.Ts = [1 - 0.009333 * x for x in self.xs]
        self.Vs = [(0.05 + 0.11 * x) * np.sqrt(1 - 0.2314 * x) for i, x in enumerate(self.xs)]
    def boundary_condition(self):
        super().boundary_condition()
        self.Ts[-1] = self.pe / self.rhos[-1]
        
if __name__ == '__main__':
    cfd_hw2 = CFD_hw2(1000)
    rhos, Vs, Ts, ts, xs = cfd_hw2.solve()
    pe = cfd_hw2.pe
    gamma = cfd_hw2.gamma
    sz = cfd_hw2.sz
    As = cfd_hw2.As
    deltaX = cfd_hw2.deltaX
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
    print(rhos)
    print(Vs)
    plt.plot(xs, rhos)
    plt.plot(xs, Vs, )
    plt.plot(xs, Ts)
    plt.plot(xs, Mas, 'o')
    plt.plot(xs, Ma_exact, 'o')
    plt.plot(xs, [Ts[i] * rhos[i] for i in range(sz)])
    plt.legend([r'$\rho$', r'$v$', r'$T$', r'Calculate $Ma$', r'Exact $Ma$', r'$c$'])
    plt.savefig('3.png')
