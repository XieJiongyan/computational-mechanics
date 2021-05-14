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
        def initial_condition(self):
            q = np.matrix([[0.5], [0], [0]])
            dq = np.matrix([[0], [0], [0]])
            return q, dq 
        def input_M(self):
            self.M = np.matrix(np.diag([1, 1, 1/12]))
            return None
        def cal_QA(self):
            g = 9.8
            return np.matrix([[0], [-g], [0]])
    k1 = K()
    k1.read_data_csv("kinematic/input/hw3_constraints.csv")
    Z1 = k1.solve()
    print(Z1.shape)

    plt.plot(ts, np.array(Z1[0, :]).flatten())
    plt.savefig('fig/CKD/x_pendulum.png')
    plt.cla()

    plt.plot(ts, np.array(Z1[3, :]).flatten())
    plt.savefig('fig/CKD/v_pendulum.png')
    plt.cla()

    plt.plot(ts, np.array(Z1[6, :]).flatten())
    plt.savefig('fig/CKD/a_pendulum.png')
    plt.cla()



