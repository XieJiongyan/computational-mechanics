import pandas as pd 
import numpy as np 
import numpy.linalg as npl 
import math 
import matplotlib.pyplot as plt 
from matplotlib import patches

class draw:
    def __init__(self):
        self.fig, self.ax = plt.subplots()
        self.r = 0.1
    def add_hinge(self, x, y):
        circle = patches.Circle([x, y], self.r, fill= True)
        self.ax.add_patch(circle)
    def xlim(self, lr):
        plt.xlim(lr)
    def ylim(self, lr):
        plt.ylim(lr)
    def cla(self):
        plt.cla()
    def add_rectangle(self, xy, width, height, angle= 0):
        xy = np.array(xy)
        angle = angle / 180 * np.pi

        xy = xy - np.array([width / 2 * math.cos(angle), width / 2 * math.sin(angle)]) + np.array([height / 2 * math.sin(angle), - height / 2 * math.cos(angle)])
        angle = angle / np.pi * 180
        rec = patches.Rectangle(xy, width, height, angle, fill= False)
        # print('angle: ', angle)
        self.ax.add_patch(rec)
    def add_pole(self, xy_start, xy_end):
        xy_start = np.array(xy_start)
        xy_end = np.array(xy_end) 
        width = npl.norm(xy_end - xy_start)
        xy = (xy_start + xy_end) / 2
        angle = np.arcsin((xy_end[1] - xy_start[1]) / width)
        angle = angle / np.pi * 180
        if xy_end[0] - xy_start[0] < 0:
            angle = 180 - angle
        # print('add_pole', xy_start[1], xy_end[1], angle)
        self.add_rectangle(xy, width, self.r, angle)
        

    def show(self):
        plt.show()
    def savefig(self, filename):
        plt.savefig(filename)


def A(phi):
    return np.matrix([[math.cos(phi), -math.sin(phi)], [math.sin(phi), math.cos(phi)]]) 
def R():
    return np.matrix([[0, -1], [1, 0]]) 
def B(phi):
    return R() * A(phi) 

def As(phii, si):
    Ai = np.matrix([[math.cos(phii), -math.sin(phii)], [math.sin(phii), math.cos(phii)]]) 
    si = si.split(' ')
    si = np.matrix([float(x) for x in si]).T 
    return Ai * si #返回列向量
def Bs(phii, si):
    return R() * As(phii, si) 

class kinematic_model:
    def __init__(self):
        self.n = 0
        self.initial_condition = 0
        self.constraints = 0 
    def read_data_from_csv(self, filename):
        df = pd.read_csv(filename)
        unique_i = pd.concat([df['i'], df['j']], ignore_index= True)
        self.constraints = df
        self.n = unique_i.value_counts().shape[0] * 3
    def read_initial_from_csv(self, filename):
        df = pd.read_csv(filename, header=None)
        if self.n != 0:
            self.n = df.shape[0]
        self.initial_condition = df.values.flatten()
    def As(self, phii, si):
        Ai = np.matrix([[math.cos(phii), -math.sin(phii)], [math.sin(phii), math.cos(phii)]]) 
        si = si.split(' ')
        si = np.matrix([float(x) for x in si]).T 
        return Ai * si #返回列向量
    def R(self):
        return np.matrix([[0, -1], [1, 0]]) 
    def Bs(self, phii, si):
        return R() * As(phii, si) 
    def cal_Phi_row(self, q, t, row):
        kind = row['kind']
        i = row['i']
        j = row['j']
        si = row['si']
        sj = row['sj']
        vi, vj = row['vi'], row['vj']
        c = row['c']

        xi = q[i * 3]
        yi = q[i * 3 + 1]
        phii = q[i * 3 + 2]
        ri = np.array([xi, yi])
        if j == j:
            j = int(j)
            xj = q[j * 3]
            yj = q[j * 3 + 1]
            rj = np.array([xj, yj])
            phij = q[j * 3 + 2]
        if kind == 'aphid':
            c = str(row['c']).replace("pi", "np.pi").replace('^', '**')
            # print('kind 105:', c, phii - eval(c), phii, q)
            return phii - eval(c)
        elif kind == 'ax':
            c = float(c)
            return xi + As(phii, si)[0, 0] - c
        elif kind == 'ay':
            c = float(c)
            return yi + As(phii, si)[1, 0] - c 
        elif kind == 'r':
            Aisi = np.array(As(phii, si))
            Ajsj = np.array(As(phij, sj))
            return (rj + Ajsj.T - ri - Aisi.T).flatten()
        elif kind == 't':
            vvi = vi.split(' ')
            vvi = np.matrix([float(x) for x in vvi]).T 
            Aisi = np.array(As(phii, si))
            Ajsj = np.array(As(phij, sj))
            h = np.matrix(rj + Ajsj.T - ri - Aisi.T).T
            # print(Bs(phii, vi).T, h)
            phi1 = (Bs(phii, vi).T * h)[0, 0]
            phi2 = (- vvi.T * Bs(phij - phii, vj))[0, 0]
            return np.array([phi1, phi2])
        else:
            print('Constrant kind unrecognized.')
            return 0
    def cal_Phi(self, q, t, constraints):
        Phi = np.array([])
        for index, row in constraints.iterrows():
            i, j = row['i'], row['j']
            # print('134:', self.cal_Phi_row(q, t, row))
            Phi = np.hstack([Phi, self.cal_Phi_row(q, t, row)])
        return np.matrix(Phi).T # 返回一个列向量
    def cal_fq_row(self, q, t, row):
        kind = row['kind']
        i = row['i']
        j = row['j']
        si = row['si']
        sj = row['sj']
        vi, vj = row['vi'], row['vj']
        c = row['c']

        xi = q[i * 3]
        yi = q[i * 3 + 1]
        phii = q[i * 3 + 2]
        ri = np.array([xi, yi])
        if j == j:
            j = int(j)
            xj = q[j * 3]
            yj = q[j * 3 + 1]
            rj = np.array([xj, yj])
            phij = q[j * 3 + 2]

        szq = len(q)
        Phi_q = 0
        if kind == 'aphid':
            Phi_q = np.zeros([1, szq])
            Phi_q[0][i * 3 + 2] = 1
            return Phi_q 
        elif kind == 'ax':
            Phi_q = np.zeros([1, szq])
            Phi_q[0, i * 3] = 1
            Phi_q[0, i * 3 + 2] = Bs(phii, si)[0, 0]
        elif kind == 'ay':
            Phi_q = np.zeros([1, szq])
            Phi_q[0, i * 3 + 1] = 1
            Phi_q[0, i * 3 + 2] = Bs(phii, si)[1, 0]
        elif kind == 'r':
            Phi_q = np.zeros([2, szq])
            Phi_q[:, [i * 3, i * 3 + 1]] = -np.identity(2)
            Phi_q[:, [j * 3, j * 3 + 1]] = np.identity(2) 
            # print(Phi_q[:, i * 3 + 2], Bs(phii, si))
            Phi_q[:, i * 3 + 2] = -Bs(phii, si).flatten()
            Phi_q[:, j * 3 + 2] = Bs(phij, sj).flatten()
        elif kind == 't':
            Phi_q = np.zeros([2, szq])
            vvi = vi.split(' ')
            vvi = np.matrix([float(x) for x in vvi]).T 
            Phi_q[0, [i * 3, i * 3 + 1]] = -(Bs(phii, vi)).T
            Phi_q[1, [i * 3, i * 3 + 1]] = np.array([[0, 0]])
            rirj = np.matrix((rj - ri).reshape([2, 1]))
            Phi_q[0, i * 3 + 2] = - (As(phii, vi)).T * rirj - vvi.T * As(phij - phii, sj)
            Phi_q[1, i * 3 + 2] = - vvi.T * As(phij - phii, vj) 
            Phi_q[0, [j * 3, j * 3 + 1]] = (Bs(phii, vi)).T 
            Phi_q[1, [j * 3, j * 3 + 1]] = np.array([[0, 0]]) 
            Phi_q[0, j * 3 + 2] = vvi.T * As(phij - phii, sj) 
            Phi_q[1, j * 3 + 2] = vvi.T * As(phij - phii, vi)
        else:
            print('''kind is wrong''')
        return np.matrix(Phi_q) 
    def cal_fq(self, q, t, constraints):
        Phi_q = np.empty((0, self.n))
        for index, row in constraints.iterrows():
            Phi_q = np.vstack([Phi_q, self.cal_fq_row(q, t, row)])
        return Phi_q
    def cal_gamma_row(self, q, dq, t, row):
        kind = row['kind']
        i = row['i']
        j = row['j']
        si = row['si']
        sj = row['sj']
        vi, vj = row['vi'], row['vj']
        c = row['c']

        xi = q[i * 3]
        yi = q[i * 3 + 1]
        phii = q[i * 3 + 2]
        dphii = dq[i * 3 + 2]
        ri = np.matrix([xi, yi]).reshape([2, 1])
        dxi = dq[i * 3]
        dyi = dq[i * 3 + 1]
        dri = np.matrix([[dxi, dyi]]).T
        if j == j:
            j = int(j)
            xj = q[j * 3]
            yj = q[j * 3 + 1]
            rj = np.matrix([xj, yj]).reshape([2, 1]) #列向量
            phij = q[j * 3 + 2]
            dphij = dq[j * 3 + 2]
            dxj = dq[j * 3]
            dyj = dq[j * 3 + 1]
            drj = np.matrix([[dxj, dyj]]).T
        if kind == 'aphid':
            c = str(row["c''"]).replace("pi", "np.pi").replace('^', '**')
            # print(eval(c))
            return eval(c)
        elif kind == 'ax':
            return (As(phii, si) * dphii ** 2)[0, 0]
        elif kind == 'ay':
            return (As(phii, si) * dphii ** 2)[1, 0]
        elif kind == 'r':
            return np.array(As(phij, sj) * dphij ** 2 - As(phii, si) * dphii ** 2).flatten()
        elif kind == 't':
            print(rj)
            t1 = (dphii ** 2 * Bs(phii, vi).T * (rj - ri) + 2 * dphii * As(phii, vi).T * (drj - dri))[0, 0]
            return np.array([t1, 0])
        else:
            return 0
    def cal_gamma(self, q, dq, t):
        gamma = np.array([])
        # print('origin', gamma)
        for index, row in self.constraints.iterrows():
            # print('1', gamma)
            gamma = np.hstack([gamma, self.cal_gamma_row(q, dq, t, row)])
            # print(gamma)
        return np.matrix(gamma).T # 返回列向量
    def cal_v_row(self, q, t, row):
        kind = row['kind']
        if kind == 'aphid':
            c = str(row["c'"]).replace("pi", "np.pi").replace('^', '**')
            # print(eval(c))
            return eval(c)
        elif kind == 'r':
            return [0, 0]
        elif kind == 't':
            return [0, 0]
        else:
            return 0
    def cal_v(self, q, t):
        v = np.array([])
        for index, row in self.constraints.iterrows():
            v = np.hstack([v, self.cal_v_row(q, t, row)])
        return np.matrix(v).T
    # @param ts: 运行时段
    # Return [phi, phi, a]
    def directly_solve(self, ts):
        if type(ts) != list:
            print("Wrong input ts, not a list")
        dt = ts[1] - ts[0]
        q = [float(x) for x in self.initial_condition]
        n = len(q)
        constraints = self.constraints
        if n != self.n:
            print('n is wrong')
        Z = np.zeros([ts.shape[0], n * 3])
        for it, t in enumerate(ts):
            print(it, q)
            Phi = self.cal_Phi(q, t, constraints)
            Phi_q = self.cal_fq(q, t, constraints)
            i = 1
            while npl.norm(self.cal_Phi(q, t, constraints)) > 1e-6:
                Phi = self.cal_Phi(q, t, constraints)
                Phi_q = self.cal_fq(q, t, constraints)
                if abs(npl.det(Phi_q)) < 1e-4:
                    print('Improper initial value')
                    print(Phi_q) 
                    print(npl.det(Phi_q)) 
                q = np.matrix(np.array(q).reshape([n, 1]))
                q -= npl.inv(Phi_q) * Phi 
                q = np.array(q).reshape(n, 1).flatten().tolist()
                # print('291:', q)
                i += 1 
                if i == 100:
                    print('Improper initial value, i >= maxi')
                    print(t)
                    return Z

            dq = npl.inv(Phi_q) * self.cal_v(q, t) 
            q = np.array(q).reshape(n, 1)
            # print(q, dq)
            # print(np.array(dq).flatten())
            gamma = self.cal_gamma(q, np.array(dq).flatten(), t)
            # print(gamma)
            ddq = npl.inv(Phi_q) * gamma
            q = np.matrix(q)
            # print('q: ', q)
            # print('dq * dt: ', dq * dt)
            # print('0.5 * dt ** 2 * ddq :', 0.5 * dt ** 2 * ddq )
            q += dq * dt + 0.5 * dt ** 2 * ddq 
            q = np.array(q).reshape(n, 1).flatten().tolist()
            # print('q:', q)
            dq = np.array(dq).flatten()
            ddq = np.array(ddq).flatten()

            Z[it, :] = np.hstack([np.array(q), dq, ddq])
        return Z

class kinematic_base:
    # 在本class中，数字数据只有两种类型：标量(float)，矩阵(np.matrix)
    def __init__(self):
        pass 
    def read_data_csv(self, filename):
        df = pd.read_csv(filename)
        unique_i = pd.concat([df['i'], df['j']], ignore_index= True)
        self.constraints = df
        self.n = unique_i.value_counts().shape[0] * 3
        self.q = np.matrix(np.zeros([self.n, 1])) 
        self.dq = np.matrix(np.zeros([self.n, 1])) 
        self.ddq = np.matrix(np.zeros([self.n, 1]))
    def get_values(self, row):
        self.kind = row['kind']

        # get c, c', c'', their types are all str.
        if self.kind == 'ax' or self.kind == 'ay':
            self.c = float(row['c'])
        else:
            if type(row['c']) == str:
                self.c = str(row["c"]).replace("pi", "np.pi").replace('^', '**')
            if type(row['c']) == str:
                self.c_ = str(row["c'"]).replace("pi", "np.pi").replace('^', '**')
            if type(row['c']) == str:
                self.c__ = str(row["c''"]).replace("pi", "np.pi").replace('^', '**')



        # get information of i
        self.i = row['i']
        if type(row['si']) == str:
            self.si = np.matrix([float(x) for x in row['si'].split(' ')]).T 
        if type(row['vi']) == str:
            self.vi = np.matrix([float(x) for x in row['vi'].split(' ')]).T 
        c = row['c']

        self.xi = self.q[self.i * 3, 0] # scalar
        self.yi = self.q[self.i * 3 + 1, 0] # scalar
        self.phii = self.q[self.i * 3 + 2, 0] # scalar
        self.ri = np.matrix([[self.xi], [self.yi]])
        self.dxi = self.dq[self.i * 3, 0] # scalar
        self.dyi = self.dq[self.i * 3 + 1, 0] # scalar 
        self.dri = np.matrix([[self.dxi], [self.dyi]])
        self.dphii = self.dq[self.i * 3 + 2, 0] # scalar

        # get information of j
        self.j = row['j']
        if not np.isnan(self.j):
            self.j = int(self.j)
            if not np.isnan(row['sj']):
                self.sj = np.matrix([float(x) for x in row['sj'].split(' ')]).T 
            if not np.isnan(row['vj']):
                self.vj = np.matrix([float(x) for x in row['vj'].split(' ')]).T 
            self.xj = self.q[self.j * 3, 0]
            self.yj = self.q[self.j * 3 + 1, 0]
            self.phij = self.q[self.j * 3 + 2, 0]
            self.rj = np.matrix([[self.xj], [self.yj]])
            self.dxj = self.dq[self.j * 3, 0]
            self.dyj = self.dq[self.j * 3 + 1, 0]
            self.drj = np.matrix([[self.dxj], [self.dyj]])
            self.dphij = self.dq[self.j * 3 + 2, 0]
    def cal_Phi_row(self, row):
        if self.kind == 'aphid':
            return self.phii - eval(self.c)
        elif self.kind == 'ax':
            return self.xi + (A(self.phii) * self.si)[0, 0] - self.c
        elif self.kind == 'ay':
            return self.yi + (A(self.phii) * self.si)[1, 0] - self.c 
        elif self.kind == 'r':
            return self.rj + A(self.phij) * self.sj - self.ri - A(self.phii) * self.si 
        elif self.kind == 't':
            h = rj + A(self.phij) * sj - ri - A(self.phii) * si 
            phi1 = (vi.T * B(self.phii).T * h)[0, 0] #scalar
            phi2 = (-vi.T * B(self.phij - self.phii) * vj)[0, 0]
            return np.matrix([[phi1], [phi2]])
        elif self.kind == 'aphi':
            return self.phii - self.c
        else:
            print('Constrant kind unrecognized.')
            return 0
    def cal_Phi(self):
        Phi = np.matrix(np.empty([0, 1]))
        for index, row in self.constraints.iterrows():
            i, j = row['i'], row['j']
            self.get_values(row)
            Phi = np.vstack([Phi, self.cal_Phi_row(row)])
        return Phi # 返回一个列向量
    def cal_fq_row(self, row):
        szq = self.q.shape[0] #scalar
        Phi_q = 0
        if self.kind == 'aphid':
            Phi_q = np.matrix(np.zeros([1, szq]))
            Phi_q[0, self.i * 3 + 2] = 1
            return Phi_q 
        elif self.kind == 'ax':
            Phi_q = np.matrix(np.zeros([1, szq]))
            Phi_q[0, self.i * 3] = 1
            Phi_q[0, self.i * 3 + 2] = (B(self.phii) * self.si)[0, 0]
        elif self.kind == 'ay':
            Phi_q = np.matrix(np.zeros([1, szq]))
            Phi_q[0, self.i * 3 + 1] = 1
            Phi_q[0, self.i * 3 + 2] = (B(self.phii) * self.si)[1, 0]
        elif self.kind == 'r':
            Phi_q = np.matrix(np.zeros([2, szq]))
            Phi_q[:, [self.i * 3, self.i * 3 + 1]] = np.matrix(-np.identity(2))
            Phi_q[:, [self.j * 3, self.j * 3 + 1]] = np.matrix(np.identity(2))
            Phi_q[:, self.i * 3 + 2] = - B(self.phii) * self.si 
            Phi_q[:, self.j * 3 + 2] = B(self.phij) * self.sj 
        elif self.kind == 't':
            Phi_q = np.matrix(np.zeros([2, szq]))
            Phi_q[0, [self.i * 3, self.i * 3 + 1]] = -self.vi.T * B(self.phii).T
            Phi_q[1, [self.i * 3, self.i * 3 + 1]] = np.matrix([[0, 0]])
            Phi_q[0, self.i * 3 + 2] = -self.vi.T * A(self.phii).T * (self.rj -self.ri) - self.vi.T * A(self.phij - self.phii) * self.sj 
            Phi_q[1, self.i * 3 + 2] = -self.vi.T * A(self.phij - self.phii) * self.vj 
            Phi_q[0, [self.j * 3, self.j * 3 + 1]] = self.vi.T * B(self.phii).T
            Phi_q[1, [self.j * 3, self.j * 3 + 1]] = np.matrix([[0, 0]]) 
            Phi_q[0, self.j * 3 + 2] = self.vi.T * A(self.phij - self.phii) * self.sj  
            Phi_q[1, self.j * 3 + 2] = self.vi.T * A(self.phij - self.phii) * self.vi 
        elif self.kind == 'aphi':
            Phi_q = np.matrix(np.zeros([1, szq]))
            Phi_q[0, self.i * 3 + 2] = 1
        else:
            print('''kind is wrong''')
        return Phi_q
    def cal_fq(self):
        Phi_q = np.matrix(np.empty((0, self.n)))
        for index, row in self.constraints.iterrows():
            self.get_values(row)
            Phi_q = np.vstack([Phi_q, self.cal_fq_row(row)])
        return Phi_q
    def cal_v_row(self, row):
        if self.kind == 'aphid':
            return eval(self.c_)
        elif self.kind == 'r':
            return np.matrix([[0], [0]])
        elif self.kind == 't':
            return np.matrix([[0], [0]])
        else:
            return 0
    def cal_v(self,):
        v = np.matrix(np.empty([0, 1]))
        for index, row in self.constraints.iterrows():
            self.get_values(row)
            v = np.vstack([v, self.cal_v_row(row)])
        return v
    def cal_gamma_row(self, row):
        if self.kind == 'aphid':
            return eval(self.c__)
        elif self.kind == 'ax':
            return (A(self.phii) * self.si * self.dphii ** 2)[0, 0]
        elif self.kind == 'ay':
            return (A(self.phii) * self.si * self.dphii ** 2)[1, 0]
        elif self.kind == 'r':
            return A(self.phij) * self.sj * self.dphij ** 2 - A(self.phii) * self.si * self.dphii ** 2
        elif self.kind == 't':
            t1 = (self.dphii ** 2 * self.vi.T * B(self.phii).T * (self.rj - self.ri) + 2 * self.dphii * self.vi.T * A(self.phii).T * (self.drj - self.dri))[0, 0] 
            return np.matrix([[t1], [0]])
        elif self.kind == 'aphi':
            return 0
        else:
            print("Unknown kind to calculate gamma")
            return 0
    def cal_gamma(self):
        gamma = np.matrix(np.empty([0, 1]))
        for index, row in self.constraints.iterrows():
            self.get_values(row)
            gamma = np.vstack([gamma, self.cal_gamma_row(row)])
        return gamma # 返回列向量

class kinematic(kinematic_base):
    def __init__(self):
        super().__init__()
    def read_initial_csv(self, filename):
        df = pd.read_csv(filename, header=None)
        if self.n != 0:
            self.n = df.shape[0]
        self.initial_condition = df.values.flatten()
    # @param ts: 运行时段
    # Return [phi, phi, a]
    def solve(self, ts):
        if type(ts) != list and type(ts) != np.ndarray:
            print("Wrong input ts, not a list")
        dt = ts[1] - ts[0]
        self.q = np.array([float(x) for x in self.initial_condition]).reshape(-1, 1)
        n = len(self.q)
        constraints = self.constraints
        if n != self.n:
            print('n is wrong')
        Z = np.matrix(np.zeros([n * 3, ts.shape[0]])) # return value
        for it, t in enumerate(ts):
            self.t = t
            Phi = self.cal_Phi()
            Phi_q = self.cal_fq()
            i = 1
            while npl.norm(self.cal_Phi()) > 1e-6:
                Phi = self.cal_Phi()
                Phi_q = self.cal_fq()
                if abs(npl.det(Phi_q)) < 1e-4:
                    print('Improper initial value')
                    print(Phi_q) 
                    print(npl.det(Phi_q)) 
                self.q -= npl.inv(Phi_q) * Phi 
                i += 1 
                if i == 100:
                    print('Improper initial value, i >= maxi')
                    print(t)
                    return Z

            self.dq = npl.inv(Phi_q) * self.cal_v() 
            gamma = self.cal_gamma()
            self.ddq = npl.inv(Phi_q) * gamma
            self.q += self.dq * dt + 0.5 * dt ** 2 * self.ddq 

            Z[:, it] = np.vstack([self.q, self.dq, self.ddq])
        return Z


class positive_dynamical_problem(kinematic_base):
    def __init__(self):
        return None 
    def input_M(self):
        self.M = np.matrix(np.diag([1, 1, 1/12]))
        return None 
    def input_t(self):
        self.deltat = 0.001 
        self.te = 5 # 结束时间
    def cal_QA(self):
        g = 9.8
        return np.matrix([[0], [-g], [0]])
    def f(self, t, y): # return dy
        self.t = t 
        self.y = y
        dy = np.matrix([[0]] * 6)
        self.q = self.y[0:self.n, 0]
        self.dq = self.y[self.n:self.n * 2, 0] 
        QA = self.cal_QA()

        fq = self.cal_fq()
        gamma = self.cal_gamma()

        lines2 = fq.shape[0]
        zeros = np.matrix(np.zeros([lines2, lines2]))
        A = np.vstack([np.hstack([self.M, fq.T]), np.hstack([fq, zeros])])
        B = np.vstack([QA, gamma])

        X = npl.inv(A) * B 
        q__, lbd = X[0 : self.n, 0], X[self.n:, 0]

        f = np.vstack([self.dq, q__])
        dy = f

        return f, lbd
    def initial_condition(self):
        # q: q, dq: q对时间的导数 都是np.matrix
        q = np.matrix([[0.5], [0], [0]])
        dq = np.matrix([[0], [0], [0]])
        return q, dq 
    def initial_y(self):
        return None 
    def step(self, t):
        f = [0] * 4
        ks = [0, 0.5, 0.5, 1]

        origin_t = self.t 
        origin_y = self.y
        for i in range(4):
            f[i], _ = self.f(origin_t + ks[i] * self.deltat, origin_y + ks[i] * self.deltat * f[i - 1]) # no bug when i = 0
        self.y = origin_y + 1/6 * (f[0] + 2 * f[1] + 2 * f[2] + f[3]) * self.deltat 
        self.t = origin_t
    def solve(self):
        self.input_M()
        self.input_t()
        self.q, self.dq = self.initial_condition()
        self.y = np.vstack([self.q, self.dq]) 
        self.status = np.matrix(np.empty([self.n * 3, 0]))
        for t in np.arange(0, self.te, self.deltat):
            self.t = t
            self.step(t)
            self.constraint_stablization() 
            dy, _ = self.f(t + self.deltat, self.y)
            self.status = np.hstack([self.status, np.vstack([self.y, dy[self.n:self.n * 2, 0]])])
        return self.status 
    def constraint_stablization(self):
        Phiq = self.cal_fq()
        Phi = self.cal_Phi()
        self.y[0:self.n, 0] -= Phiq.T * npl.inv(Phiq * Phiq.T) * Phi
        
