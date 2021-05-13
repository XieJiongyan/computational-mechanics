import kinematics
import math
import numpy as np 
import matplotlib.pyplot as plt 


if __name__ == '__main__':
    deltat = 0.01
    te = 1.5
    ts = np.arange(0, te, deltat)

    class K(kinematics.positive_dynamical_problem):
        def input_M(self):
            self.M = np.matrix(np.diag([200, 200, 450, 35, 35, 35, 25, 25, 0.02]))
        def input_t(self):
            self.deltat = deltat 
            self.te = te # 结束时间
        def cal_QA(self):
            g = 9.8
            x3 = self.q[2 * 3, 0]
            dx3 = self.q[2 * 3 + 1, 0] 
            if dx3 < 0:
                F = 0 
            elif x3 < 1.5:
                F = 0
            elif x3 < 5:
                F = - 282857 / (6 - x3) + 62857 
            else:
                F = -110000 * (1 - math.sin(2 * 3.1415926 * (x3 - 5.25)))
            return np.matrix([[0], [-200 * g], [41450], [0], [-35 * g], [0], [F], [-25 * g], [0]])
        def initial_condition(self):
            # q: q, dq: q对时间的导数 都是np.matrix
            q = np.matrix([[0], [1], [np.pi], [1.435], [1], [-0.608], [2.47], [0], [0]])
            dq = np.matrix([[-30], [0], [30], [-60], [0], [0], [-60], [0], [0]])
            return q, dq 

    k1 = K()
    k1.read_data_csv("kinematic/input/Bighw2_1_constraints.csv")
    Z1 = k1.solve()
    print(Z1[0, :])

    plt.plot(ts, np.array(Z1[2, :]).flatten())
    plt.savefig('fig/CKD/Bighw2_1_phi.png')
    plt.cla()

    plt.plot(ts, np.array(Z1[9 + 2, :]).flatten())
    plt.savefig('fig/CKD/Bighw2_1_omega.png')
    plt.cla()

