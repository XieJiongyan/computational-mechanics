import CSM 

if __name__ == '__main__':
    csm3_7 = CSM.CSM_trass()
    csm3_7.read('CSM/book3-7.txt')
    csm3_7.solve() 
    ds, inside_F = csm3_7.get_output()
    print(ds)
    for i, F in enumerate(inside_F):
        print("F" + str(i) + " = ")
        print(F)