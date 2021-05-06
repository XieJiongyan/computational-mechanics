if __name__ == '__main__':
    import kinematics
    import math
    import numpy as np 
    import matplotlib.pyplot as plt 

    hw3 = kinematics.positive_dynamical_problem() 
    hw3.read_data_from_csv("kinematic/hw3_constraints.csv")
    status = hw3.solve() 
    print(status[-1, :]) 
    print(status[:, 0])
    sz = status.shape[0] 
    plt.plot(range(sz), status[:, 0])
    plt.plot(range(sz), status[:, 1])
    plt.plot(range(sz), status[:, 2])
    plt.legend(['x', 'y', 'phi'])
    plt.savefig('fig/CKD/hw3_1.png')

    plt.clf()
    plt.plot(range(sz), status[:, 3])
    plt.plot(range(sz), status[:, 4])
    plt.plot(range(sz), status[:, 5])
    plt.legend(['dx', 'dy', 'dphi'])
    plt.savefig('fig/CKD/hw3_2.png')

    plt.clf()
    plt.plot(range(sz), status[:, 6])
    plt.plot(range(sz), status[:, 7])
    plt.plot(range(sz), status[:, 8])
    plt.legend(['ddx', 'ddy', 'ddphi'])
    plt.savefig('fig/CKD/hw3_3.png')
