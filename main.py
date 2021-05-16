import math
import numpy as np 
from time import time


def compress(x,flip = False):
    #linearly independent combination of differences
    if flip == False:
        a1 =  2*(x[1] - x[0])
        a2 =  2*(x[2] - x[1])
        a3 =  2*(x[3] - x[2])
    elif flip == True:
        a1 =  2*(x[0] - x[1])
        a2 =  2*(x[1] - x[2])
        a3 =  2*(x[2] - x[3])

    return np.array([a1,a2,a3])

def square_sum_idx(x,y= None,z = None,idx = 0, num=3):
    if num == 3:
        return x[idx]**2 + y[idx]**2 + z[idx]**2
    elif num == 2:
        return x[idx]**2 - x[idx+1]**2 
    else: 
        return x[idx]**2

def get_e(x,y,z,p):
    e1 =  square_sum_idx(x,y,z,0) - square_sum_idx(x,y,z,1) - square_sum_idx(p,idx = 0, num =2) 
    e2 =  square_sum_idx(x,y,z,1) - square_sum_idx(x,y,z,2) - square_sum_idx(p,idx = 1, num =2)
    e3 =  square_sum_idx(x,y,z,2) - square_sum_idx(x,y,z,3) - square_sum_idx(p,idx = 2, num =2)
    return np.array([e1,e2,e3])

def make_matrix(a,b,c):

    M = np.dstack((a,b,c)).reshape((3,3))
    return  M

def get_a(MI, MG,x,y,z,pt, i = 0):
    #final compression
    def sub_alp1(MI, MG, x, idx):
        return MI[idx]*(MG[idx] - x[i])
    def sub_alp0(MG,x,idx):
        return (MG[idx] - x[i])**2

    a0 = sub_alp0(MG,x,0) + sub_alp0(MG,y,1) + sub_alp0(MG,z,2) - pt[i]**2
    a1 = 2*(pt[i] + sub_alp1(MI, MG, x, 0) + sub_alp1(MI, MG, y, 1) + sub_alp1(MI, MG, z, 2))
    a2 = (MI[0]**2) + (MI[1]**2) + (MI[2]**2) - 1

    return np.array([a0,a1,a2])

def get_b(alpha):
    temp = (alpha[1]**2 - 4*alpha[2]*alpha[0])
    b = (-alpha[1] + math.sqrt(temp)) / (2*alpha[2])
    return b

def save_vars(x,y,z,p):
    np.savetxt('x.txt', x)
    np.savetxt('y.txt', y)
    np.savetxt('z.txt', z)
    np.savetxt('p.txt', p)

def load_vars():
    x=np.loadtxt('x.txt')
    y=np.loadtxt('y.txt')
    z=np.loadtxt('z.txt')
    p=np.loadtxt('p.txt')
    return x,y,z,p

def four_gps_solwithsp(x,y,z,pt):
    "4 satellite spherical-plane algorithm"
    alpha = compress(x)
    beta = compress(y)
    gamma = compress(z)
    delta = compress(pt,True)
    epsilon =  get_e(x,y,z,pt)
    
    #Define Matrixes
    H = make_matrix(alpha,beta,gamma)
    
    #Matrix inverse for division
    H_inv = np.linalg.pinv(H)
    
    #vectors (from 5)
    B = 1*H_inv@(delta.T) 
    C = -1*H_inv@(epsilon.T)
    
    b_i = []
    x_list = []
    y_list = []
    z_list = []
    for i in range(4):
        a = get_a(B,C,x,y,z,pt,i)
        b = get_b(a)
        x_i, y_i, z_i =B*b + C
        b_i.append(b)
        x_list.append(x_i)
        y_list.append(y_i)
        z_list.append(z_i)
    b=np.array(b_i)
    x_i=np.array(x_list)
    y_i=np.array(y_list)
    z_i=np.array(z_list)
    
    # print('a2,a1,a1:')
    # print(np.flip(a,0))
    # print('b:')
    # print(b)
    # print('x,y,z:')
    # print(x_i, y_i, z_i)
    return x_i, y_i, z_i, b

def _verification(x,x_i, y,y_i, z,z_i,p_i, b,n):
    p_hat = []
    for i in range(n):
        rho_i = math.sqrt((x-x_i[i])**2 + (y-y_i[i])**2 + (z-z_i[i])**2)
        p = rho_i + b
        p_hat.append(p)
    p_hat = np.array(p_hat)
    score = np.linalg.norm((p_i-p_hat))
    print(score)
def verification(x,x_i, y,y_i, z,z_i,p_i, b,n):
    p_hat = []
    for i in range(n):
        rho_i = math.sqrt((x[i]-x_i[i])**2 + (y[i]-y_i[i])**2 + (z[i]-z_i[i])**2)
        p = rho_i + b[i]
        p_hat.append(p)
    p_hat = np.array(p_hat)
    error = np.linalg.norm((p_i-p_hat))
    if error < 10**(-6):
        print(error)
        print('Reasonable Error')



def main():
    #Define GPS values
    tic = time()
    x_i,y_i,z_i,p_i = load_vars()
    n = x_i.shape[0]
    # _four_gps_sol(x,y,z,pt)
    x,y,z,b = four_gps_solwithsp(x_i,y_i,z_i,p_i)

    verification(x,x_i, y,y_i, z,z_i,p_i, b,n)
    toc = time()
    print(toc -tic)  

if __name__ ==  "__main__":
    main()



