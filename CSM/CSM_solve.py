import CSM 
import numpy as np 

if __name__ == '__main__':
    # csm3_7 = CSM.CSM_trass('CSM/book3-7.txt')
    # ds, inside_F, out_F = csm.get_output()
    # print(ds)
    # for i, F in enumerate(inside_F):
    #     print("F" + str(i) + " = ")
    #     print(F)

    csm = CSM.CSM_trass('CSM/tao.txt')
    ds, inside_F, out_F = csm.get_output()
    print(np.round(csm.K, 3))
    print(ds)
    for i, F in enumerate(inside_F):
        print("F" + str(i) + " = ")
        print(F)
    for i, F in enumerate(out_F):
        if i % 2 == 0:
            print("Fx" + str(i//2) + " = ")
            print(F)
        else:
            print("Fy" + str(i//2) + " = ")
            print(F)
