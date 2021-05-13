import CFD 
import numpy as np
import numpy.linalg as npl
import matplotlib.pyplot as plt
import math

if __name__ == '__main__':
    cfd_hw1 = CFD.CFD_Mac_conservation(1)
    rhos, Vs, Ts, ts, xs = cfd_hw1.solve()
    print(len(rhos))
    print(rhos)
    print(Vs)
    plt.plot(xs, rhos)
    plt.savefig('fig/CFD/hw3.png')
