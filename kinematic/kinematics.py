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
    def directly_solve(self):
        t0 = 0
        te = 2
        num = 50
        ts = np.linspace(t0, te, num + 1)
        q = [float(x) for x in self.initial_condition]
        n = len(q)
        constraints = self.constraints
        if n != self.n:
            print('n is wrong')
        Z = np.zeros([num + 1, n])
        for it, t in enumerate(ts):
            i = 1
            while npl.norm(self.cal_Phi(q, t, constraints)) > 1e-6:
                Phi = self.cal_Phi(q, t, constraints)
                Phi_q = self.cal_fq(q, t, constraints)
                if abs(npl.det(Phi_q)) < 1e-4:
                    print('Improper initial value')
                    print(Phi_q) 
                    print(npl.det(Phi_q)) 
                dq = - npl.inv(Phi_q) * Phi 
                q = np.array(q).reshape(n, 1)
                # print(q, dq)
                q += dq
                q = q.flatten().tolist()
                i += 1 
                # print('times: ', t, i, self.cal_Phi(q, t, constraints), Phi_q, dq, q)
                if i == 10:
                    print('Improper initial value, i >= maxi')
                    print(t)
                    return Z
            Z[it, :] = np.array(q)
        return Z

n = 0
initial_condition = 0
constraints = 0

def read_data_from_csv(filename):
    df = pd.read_csv(filename)
    # print(df.at[0, 'kind'])
    unique_i = pd.concat([df['i'], df['j']], ignore_index= True)
    global n, constraints
    constraints = df
    n = unique_i.value_counts().shape[0] * 3
    # print(n)
    # print(unique_i)

def read_initial_from_csv(filename):
    df = pd.read_csv(filename, header=None)
    global n, initial_condition
    if n != 0:
        n = df.shape[0]
    initial_condition = df.values.flatten()
    # print(initial_condition)

def As(phii, si):
    Ai = np.matrix([[math.cos(phii), -math.sin(phii)], [math.sin(phii), math.cos(phii)]]) 
    si = si.split(' ')
    si = np.matrix([float(x) for x in si]).T 
    return Ai * si #返回列向量

def R():
    return np.matrix([[0, -1], [1, 0]]) 

def Bs(phii, si):
    return R() * As(phii, si) 

def cal_Phi_row(q, t, row):
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
        return 0

def cal_Phi(q, t, constraints):
    Phi = np.array([])
    for index, row in constraints.iterrows():
        i, j = row['i'], row['j']
        Phi = np.hstack([Phi, cal_Phi_row(q, t, row)])
    return np.matrix(Phi).T # 返回一个列向量

def cal_fq_row(q, t, row):
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
        '''kind is wrong'''
    return np.matrix(Phi_q) 

def cal_fq(q, t, constraints):
    global n
    Phi_q = np.empty((0, n))
    for index, row in constraints.iterrows():
        Phi_q = np.vstack([Phi_q, cal_fq_row(q, t, row)])
    return Phi_q
    
def directly_solve():
    t0 = 0
    te = 2
    num = 50
    ts = np.linspace(t0, te, num + 1)
    q = [float(x) for x in initial_condition]
    n = len(q)
    Z = np.zeros([num + 1, n])
    for it, t in enumerate(ts):
        i = 1
        while npl.norm(cal_Phi(q, t, constraints)) > 1e-6:
            Phi = cal_Phi(q, t, constraints)
            Phi_q = cal_fq(q, t, constraints)
            if abs(npl.det(Phi_q)) < 1e-4:
                print('Improper initial value')
                print(Phi_q) 
                print(npl.det(Phi_q)) 
            dq = - npl.inv(Phi_q) * Phi 
            q = np.array(q).reshape(n, 1)
            # print(q, dq)
            q += dq
            q = q.flatten().tolist()
            i += 1 
            # print('times: ', t, i, cal_Phi(q, t, constraints), Phi_q, dq, q)
            if i == 10:
                print('Improper initial value, i >= maxi')
                print(t)
                return Z
        Z[it, :] = np.array(q)
    return Z

