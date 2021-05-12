import CSM 

if __name__ == '__main__':
    # csm3_7 = CSM.CSM_trass()
    # csm3_7.read('CSM/book3-7.txt')
    # csm3_7.solve() 
    # ds, inside_F = csm3_7.get_output()
    # print(ds)
    # for i, F in enumerate(inside_F):
    #     print("F" + str(i) + " = ")
    #     print(F)

    csm = CSM.CSM_trass()
    csm.read('CSM/test.txt')
    csm.solve() 
    ds, inside_F, out_F = csm.get_output()
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
