import kinematics
import math
import numpy as np 
import matplotlib.pyplot as plt 

if __name__ == '__main__':

    t0 = 0
    te = 2
    num = 50
    ts = np.linspace(t0, te, num + 1)
    k1 = kinematics.kinematic()
    k1.read_data_csv("kinematic/input/hw1_constraints.csv")
    k1.read_initial_csv("kinematic/input/hw1_IC.csv")
    Z1 = k1.solve(ts)
    print(Z1)

    plt.plot(ts, np.array(Z1[0, :]).flatten().tolist())
    plt.savefig('fig/CKD/x1_new.png')
    plt.cla()


