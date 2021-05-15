import math
import numpy as np 


def compress(x,flip = False):
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

def get_alpha(MI, MG,x,y,z,pt):
    #final compression
    def sub_alp1(MI, MG, x, idx):
        return MI[idx]*(MG[idx] - x[3])
    def sub_alp0(MG,x,idx):
        return (MG[idx] - x[3])**2

    alpha0 = sub_alp0(MG,x,0) + sub_alp0(MG,y,1) + sub_alp0(MG,z,2) - pt[3]**2
    alpha1 = 2*(pt[3] + sub_alp1(MI, MG, x, 0) + sub_alp1(MI, MG, y, 1) + sub_alp1(MI, MG, z, 2))
    alpha2 = (MI[0]**2) + (MI[1]**2) + (MI[2]**2) - 1

    return np.array([alpha0,alpha1,alpha2])

def get_Rc(alpha):
    temp = (alpha[1]**2 - 4*alpha[2]*alpha[0])
    Rc = (-alpha[1] + math.sqrt(temp)) / (2*alpha[2])
    return Rc

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

def main():
    #Define GPS values
    x,y,z,pt = load_vars()

    a = compress(x)
    b = compress(y)
    c = compress(z)
    d = compress(pt,True)
    e =  get_e(x,y,z,pt)
    
    #Define Matrixes
    M = make_matrix(a,b,c)
    
    #Matrix inverse for division
    M_inv = np.linalg.pinv(M)
    
    #vectors
    MI = -1*M_inv@(d.T) 
    MG = -1*M_inv@(e.T)

    alpha = get_alpha(MI,MG,x,y,z,pt)
    Rc = get_Rc(alpha)
    print(MI*Rc + MG)
    print(Rc)

if __name__ ==  "__main__":
    main()



