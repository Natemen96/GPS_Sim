import math
import numpy as np 
from time import time


def compress(x,n=4,flip = False):
    "linearly independent combination of differences, equations 3"
    a_list = []
    if flip == False:
        for j in range(n-1):
            a_j = 2*(x[j+1] - x[j])
            a_list.append(a_j)
    elif flip == True:
        for j in range(n-1):
            a_j = 2*(x[j] - x[j+1])
            a_list.append(a_j)
    return np.array(a_list)

def square_sum_idx(x,y= None,z = None,j = 0, num=3):
    "equations 3"
    if num == 3:
        return x[j]**2 + y[j]**2 + z[j]**2
    elif num == 2:
        return x[j]**2 - x[j+1]**2 
    else: 
        return x[j]**2

def get_e(x,y,z,p,n=4):
    e_list = []
    for j in range(n-1):
        e_j = square_sum_idx(x,y,z,j) - square_sum_idx(x,y,z,j+1) - square_sum_idx(p,j = j, num =2) 
        e_list.append(e_j)
    return np.array(e_list)

def make_matrix(M):

    H = np.dstack(M)
    H = H.reshape(H.shape[1],H.shape[2])
    return  H

def get_a(B, C,x,y,z,pt, i = 3):
    #final compression
    def sub_alp1(B, C, x, idx):
        return B[idx]*(C[idx] - x[i])
    def sub_alp0(C,x,idx):
        return (C[idx] - x[i])**2

    a2 = (B[0]**2) + (B[1]**2) + (B[2]**2) - 1
    a1 = 2*(pt[i] + sub_alp1(B, C, x, 0) + sub_alp1(B, C, y, 1) + sub_alp1(B, C, z, 2))
    a0 = sub_alp0(C,x,0) + sub_alp0(C,y,1) + sub_alp0(C,z,2) - pt[i]**2

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

def load_vars(n = 8):
    assert n >= 4, 'n must be <= 4'
    x=np.loadtxt('x.txt')[:n]
    y=np.loadtxt('y.txt')[:n]
    z=np.loadtxt('z.txt')[:n]
    p=np.loadtxt('p.txt')[:n]
    return x,y,z,p

def four_gps_solwithsp(x,y,z,pt):
    "4 satellite spherical-plane algorithm"
    alpha = compress(x)
    beta = compress(y)
    gamma = compress(z)
    delta = compress(pt,flip = True)
    epsilon =  get_e(x,y,z,pt)
    
    #Define Matrixes
    H = make_matrix((alpha,beta,gamma))
    
    #Matrix inverse for division
    H_inv = np.linalg.pinv(H)
    
    #vectors (from 5)

    B = -1*H_inv@(delta.T) 
    C = -1*H_inv@(epsilon.T)
    
    a = get_a(B,C,x,y,z,pt)
    print(np.flip(a,0))
    b = get_b(a)
    x_i, y_i, z_i =B*b + C
    return x_i, y_i, z_i, b


def n_gps_solwithsp(x,y,z,pt,n):
    alpha = compress(x,n)
    beta = compress(y,n)
    gamma = compress(z,n)
    delta = compress(pt,n,flip = True)
    u =  get_e(x,y,z,pt,n)
    H = make_matrix((alpha,beta,gamma,delta))

    HTH_inv = np.linalg.pinv(H.T@H)
    #I believe this is a mistake in the paper, it should be -H^Tu not H^Tu
    x_i,y_i,z_i,b =  HTH_inv@-(H.T@u) 
    return x_i,y_i,z_i,b

def verification(x,x_i, y,y_i, z,z_i,p_i, b,n):
    p_hat = []
    for i in range(n):
        # rho_i = math.sqrt((x[i]-x_i[i])**2 + (y[i]-y_i[i])**2 + (z[i]-z_i[i])**2)
        # p = rho_i + b[i]
        rho_i = math.sqrt((x-x_i[i])**2 + (y-y_i[i])**2 + (z-z_i[i])**2)
        p = rho_i + b
        p_hat.append(p)
    p_hat = np.array(p_hat)
    error = np.linalg.norm((p_i-p_hat))
    print('error:')
    print(error)
    if error < 10**(-4):
        print('Reasonable Error')
    else:
        print('Could be better')



def main():
    #Define GPS values
    tic = time()
    # x_i,y_i,z_i,p_i = load_vars(4)
    # n = x_i.shape[0]
    # x,y,z,b = four_gps_solwithsp(x_i,y_i,z_i,p_i)
    x_i,y_i,z_i,p_i = load_vars()
    n = x_i.shape[0]
    x,y,z,b = n_gps_solwithsp(x_i,y_i,z_i,p_i,n)
    print(x,y,z,b)
    verification(x,x_i, y,y_i, z,z_i,p_i, b,n)
    toc = time()
    print('time:')
    print(toc -tic)  

if __name__ ==  "__main__":
    main()



