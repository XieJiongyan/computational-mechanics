import numpy as np 
import numpy.linalg as npl 
import matplotlib.pyplot as plt 
import math 

class CSM_base:
    # 前处理
    def read(self, file):
        self.f = open(file)
        pass
    # 有限元分析主体部分
    def solve(self):
        pass 
    # 后处理
    def get_output(self):
        pass
    
class CSM_trass(CSM_base):
    def r(self):
        return int(self.f.readline().strip())
    def r_line(self, f= int):
        return np.array([f(x) for x in self.f.readline().strip().split(' ')]) 
    def read(self, file):
        super().read(file)
        self.ne = self.r() 
        self.nn = self.r() 
        self.nf = self.r()
        self.es = [] # information of start node and end node, and EA of a edge unit
        self.ns = [] # information of places of nodes 
        for i in range(self.ne): # read all unit(line unit)
            self.es += [self.r_line()] 
        for i in range(self.nn): # read all nodes 
            self.ns += [self.r_line(float)] 
        
        self.nF = self.r() # how many out Force
        self.Fs = [] 
        for i in range(self.nF):
            self.Fs += [self.r_line()] # node, direction(start from 0), value 
        self.nC = self.r() # how may constraints 
        self.Cs = [] 
        for i in range(self.nC):
            self.Cs += [self.r_line()] # node, direction(start from 0)

    def solve(self):
        self.K = np.matrix(np.zeros([self.nn * self.nf, self.nn * self.nf]))
        self.base_dir = [] 
        for i in range(self.nf):
            dir = np.matrix(np.zeros([self.nf, 1])) 
            dir[i, 0] = 1 
            self.base_dir += [dir]
        for e in self.es:
            ids = np.empty([0, 1]) 
            for i in range(self.nf):
                ids = np.vstack([ids, int(e[0] * self.nf + i)]) 
            for i in range(self.nf):
                ids = np.vstack([ids, int(e[1] * self.nf + i)]) 
            ids = ids.astype(int)
            
            dir = np.matrix((self.ns[e[1]] - self.ns[e[0]]).reshape([-1, 1]))
            Ke = e[2] / npl.norm(dir) * np.matrix([[1, -1], [-1, 1]]) 
            dir = dir / npl.norm(dir)
            lbd = np.matrix(np.empty([1, 0])) 
            for f in range(self.nf):
                lbd = np.hstack([lbd, dir.T * self.base_dir[f]]) 
            lbd_zero = np.matrix(np.zeros(lbd.shape))
            T = np.vstack([np.hstack([lbd, lbd_zero]),np.hstack([lbd_zero, lbd])])

            # print(ids)
            # print(self.K)
            # print(T)
            self.K[ids, ids.T] += T.T * Ke * T 
            # print("add")
            # print(T.T * Ke * T)
        self.F = np.matrix(np.zeros([self.nn * self.nf, 1])) 
        for F in self.Fs:
            self.F[F[0] * self.nf + F[1], 0] += F[2]
        free_nodes = np.array([True] * self.nn * self.nf) 
        free_Ks = np.array([[True] * self.nn * self.nf] * self.nn * self.nf)
        for C in self.Cs: 
            free_nodes[C[0] * self.nf + C[1]] = False 
            for i in range(self.nn * self.nf):
                free_Ks[i, C[0] * self.nf + C[1]] = False 
                free_Ks[C[0] * self.nf + C[1], i] = False
        K = self.K[free_Ks].reshape(self.nf * self.nn - self.nC, self.nf * self.nn - self.nC)
        F = self.F[free_nodes] 
        # print(self.K)
        # print(K) 
        # print(F) 
        d_ = npl.inv(K) * F # 位移
        id = 0 
        self.d = np.matrix(np.empty([0, 1])) 
        for i in free_nodes:
            if i:
                self.d = np.vstack([self.d, d_[id, 0]])
                id += 1 
            else:
                self.d = np.vstack([self.d, 0]) 
    def get_output(self):
        inside_F = [] # 杆内力
        for e in self.es:
            dir = np.matrix((self.ns[e[1]] - self.ns[e[0]]).reshape([-1, 1]))
            Ke = e[2] / npl.norm(dir) * np.matrix([[1, -1], [-1, 1]]) 
            dir = dir / npl.norm(dir)
            lbd = np.matrix(np.empty([1, 0])) 
            for f in range(self.nf):
                lbd = np.hstack([lbd, dir.T * self.base_dir[f]]) 
            lbd_zero = np.matrix(np.zeros(lbd.shape))
            T = np.vstack([np.hstack([lbd, lbd_zero]),np.hstack([lbd_zero, lbd])])
            d = np.matrix(np.empty([0, 1])) 
            for i in range(self.nf):
                d = np.vstack([d, self.d[e[0] * self.nf + i, 0]]) 
            for i in range(self.nf):
                d = np.vstack([d, self.d[e[1] * self.nf + i, 0]]) 
            inside_F += [Ke * T * d]
        return self.d, inside_F, self.K * self.d