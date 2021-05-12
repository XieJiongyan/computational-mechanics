import kinematics
import math
import numpy as np 
import matplotlib.pyplot as plt 


if __name__ == '__main__':
    deltat = 0.005
    te = 5
    ts = np.arange(0, te, deltat)

    class K(kinematics.positive_dynamical_problem):
        def input_t(self):
            self.deltat = deltat 
            self.te = te # 结束时间

    k1 = K()
    k1.read_data_csv("kinematic/input/hw3_constraints.csv")
    Z1 = k1.solve()
    # print(Z1)
    print(np.max(np.array(Z1[0, :]) ** 2 + np.array(Z1[0, :]) ** 2 - 0.5))
    print(Z1[6:9, :15])
    # print(Z1[:, [49, 50, 51]])

    plt.plot(ts, np.array(Z1[0, :]).flatten())
    plt.savefig('fig/CKD/x1_new.png')
    plt.cla()

    plt.plot(ts, np.array(Z1[3, :]).flatten())
    plt.savefig('fig/CKD/v_pendulum.png')
    plt.cla()

    plt.plot(ts, np.array(Z1[6, :]).flatten())
    plt.savefig('fig/CKD/a_pendulum.png')
    plt.cla()



