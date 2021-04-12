if __name__ == '__main__':
    import kinematics
    import math
    k = kinematics.kinematic_model()
    k.read_data_from_csv("kinematic/hw2_constraints.csv")
    k.read_initial_from_csv("kinematic/hw2_IC.csv")
    Z = k.directly_solve()
    print(Z)
    
    draw = kinematics.draw()
    # draw.add_rectangle([1, 1], 5, 5)

    NOT_PLOT_IMAGE = True
    if NOT_PLOT_IMAGE:
        exit()
    import imageio
    gif_images = []
    sz = Z.shape[0]
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
        phi0 = Z[i, 2]
        phi1 = Z[i, 5]
        phi3 = Z[i, 11]
        x3 = Z[i, 9]
        draw.add_pole([0, 0], [Z[i, 6], Z[i, 7]])
        draw.add_pole([0, -2], [0 + 4 * math.cos(phi1), -2 + 4 * math.sin(phi1)])
        draw.add_pole([x3, 2], [x3 + 1.9 * math.cos(phi3), 2 + 1.9 * math.sin(phi3)])
        draw.add_hinge(0, 0)
        draw.add_hinge(0, -2) 
        draw.add_hinge(Z[i, 6], Z[i, 7])
        draw.add_hinge(0 + 4 * math.cos(phi1), -2 + 4 * math.sin(phi1))
        draw.xlim([-4, 3])
        draw.ylim([-4, 3])
        draw.savefig('kinematic/test.png')
        gif_images.append(imageio.imread('kinematic/test.png'))
        draw.cla()
    imageio.mimsave('kinematic/test2.gif', gif_images, fps= 15)
        

