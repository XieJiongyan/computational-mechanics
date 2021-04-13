if __name__ == '__main__':
    import kinematics
    import math
    import numpy as np 
    import matplotlib.pyplot as plt 

    t0 = 0
    te = 2
    num = 50
    ts = np.linspace(t0, te, num + 1)
    k1 = kinematics.kinematic_model()
    k1.read_data_from_csv("kinematic/hw1_constraints.csv")
    k1.read_initial_from_csv("kinematic/hw1_IC.csv")
    Z1 = k1.directly_solve(ts)
    print(Z1)

    plt.plot(ts, Z1[:, 0])
    plt.savefig('x1.png')
    plt.cla()
    plt.plot(ts, Z1[:, 3])
    plt.savefig('v1.png')
    plt.cla()
    plt.plot(ts, Z1[:, 6])
    plt.savefig('a1.png')
    plt.cla()

    t0 = 0
    te = 6
    num = 120
    ts = np.linspace(t0, te, num + 1)

    k2 = kinematics.kinematic_model()
    k2.read_data_from_csv("kinematic/hw2_constraints.csv")
    k2.read_initial_from_csv("kinematic/hw2_IC.csv")
    Z2 = k2.directly_solve(ts)
    print(Z2)
    
    plt.plot(ts, Z2[:, 9])
    plt.savefig('x5.png')
    plt.cla()
    plt.plot(ts, Z2[:, 9 + 12])
    plt.savefig('v5.png')
    plt.cla()
    plt.plot(ts, Z2[:, 9 + 24])
    plt.savefig('a5.png')
    plt.cla()

    draw = kinematics.draw()
    # draw.add_rectangle([1, 1], 5, 5)

    NOT_PLOT_IMAGE = True
    if NOT_PLOT_IMAGE:
        exit()
    import imageio
    gif_images = []
    sz = Z2.shape[0]
    # for i in range(sz):
    #     phii = Z[i, 2]
    #     draw.add_pole([0, 0], [1 * math.cos(phii), 1 * math.sin(phii)])
    #     draw.add_hinge(0, 0)
    #     # draw.show()
    #     draw.xlim([-2, 2])
    #     draw.ylim([-2, 2])
    #     draw.savefig('test.png')
    #     gif_images.append(imageio.imread('test.png'))
    #     draw.cla()
    for i in range(sz):
        print('image: ', i)
        phi0 = Z2[i, 2]
        phi1 = Z2[i, 5]
        phi3 = Z2[i, 11]
        x3 = Z2[i, 9]
        draw.add_pole([0, 0], [Z2[i, 6], Z2[i, 7]])
        draw.add_pole([0, -2], [x3 + 1.9 * math.cos(phi3), 2 + 1.9 * math.sin(phi3)])
        draw.add_pole([x3, 2], [x3 + 1.9 * math.cos(phi3), 2 + 1.9 * math.sin(phi3)])
        draw.add_hinge(0, 0)
        draw.add_hinge(0, -2) 
        draw.add_hinge(Z2[i, 6], Z2[i, 7])
        draw.add_hinge(0 + 4 * math.cos(phi1), -2 + 4 * math.sin(phi1))
        draw.xlim([-4, 3])
        draw.ylim([-4, 3])
        draw.savefig('kinematic/test.png')
        gif_images.append(imageio.imread('kinematic/test.png'))
        draw.cla()
    imageio.mimsave('kinematic/test2-.gif', gif_images, fps= 15)
        

