import numpy as np


def compute_crossconj(x,y,N,Lp):
    b=x*np.conj(np.array([y]).T) #diad    
    bandmat=np.eye(N, k=0)
    for k in range(1,Lp+1):
        bandmat+=np.eye(N, k=-2*k)+np.eye(N, k=2*k)        
    c=(bandmat*b)[:,::-1]

    arr=[]
    if np.mod(N,2)==1:
        for j in range(1,int((N-1)/2+1))[::-1]:
            arr.append(np.trace(c[:-2*j,2*j:]))
        arr.append(np.trace(c))
        for j in range(1,int((N-1)/2+1)):
            arr.append(np.trace(c[2*j:,:-2*j]))
    else:
        for j in range(0,int(N/2))[::-1]:
            arr.append(np.trace(c[:-2*j-1,2*j+1:]))
        for j in range(0,int(N/2))[::-1]:
            arr.append(np.trace(c[-2*j-1:,:2*j+1]))

    arr=np.array(arr)
    return arr
            
