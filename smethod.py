import numpy as np
import scipy 
import stft 
def compute_crossconj_mat(x,y,N,Lp):
    b=x*(np.conj(np.array([y]).T)) #diad    
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

    if False:
        fig=plt.figure()
        ax=fig.add_subplot(121)
        ax.imshow(np.abs(c))
        ax=fig.add_subplot(122)
        ax.imshow(np.abs(bandmat))
        
        plt.show()
            
    arr=np.array(arr)
    return arr

def compute_crossconj(x,y,N,Lp):

    arr=[]
    for k in range(0,N):
        arr.append(np.sum(np.real(x[k-Lp:k+Lp+1]*np.conj(y[k-Lp:k+Lp+1][::-1]))))
    arr=np.array(arr)
    return arr


def tfrsm(x,y=None,Lp=6,f=None,nwindow=4,silent=0,fps=None,fpe=None,itc=None):
    import stft
    print("Under debugging. Do not use.")
    if y is None:
        stftx=stft.tfrstft(x,nwindow=nwindow)
        stfty=stftx

    if itc is None:
        itc=range(0,len(x))
  
    if f is None: 
        tfrstftx=stft.tfrstft(x,itc=itc,nwindow=nwindow)
        if y is None:
            tfrstfty=tfrstftx
        else:
            tfrstfty=stft.tfrstft(y,itc=itc,nwindow=nwindow)
        
    else:
        tfrstftx=stft.tfrstft(x,f=f,itc=itc,nwindow=nwindow)            
        if y is None:
            tfrstfty=tfrstftx
        else:
            tfrstfty=stft.tfrstft(y,f=f,itc=itc,nwindow=nwindow)


    nsamplef,nsamplet=np.shape(tfrstftx)

    if fps is None or fpe is None: 
        ks=0
        ke=nsamplef-1
        sm=np.zeros((nsamplef,nsamplet),dtype=np.complex)
    else:
        ks=np.int(np.max([1,2*fps-Lp]))-1
        ke=np.int(np.min([nsamplef,2*fpe+Lp]))
        sm=np.zeros((ke-ks+1,nsamplet),dtype=np.complex)

        #    for k in range(ks,ke):
        #        kq=k-ks
    for j in range(0,nsamplet):
#    for j in range(128,129):
        sm[:,j]=compute_crossconj(tfrstftx[:,j],tfrstfty[:,j],nsamplef,Lp)
    return sm

if __name__ == "__main__":
    import sampledata as sd
    import matplotlib.pyplot as plt

    nsamp=512
    t,x=sd.genmultifm622(nsamp)
    tfr=tfrsm(x,Lp=12,nwindow=8)
    tfrstft=stft.tfrstft(x,nwindow=8)
    
    print(np.max(tfr),np.min(tfr))
    fig=plt.figure()
    ax=fig.add_subplot(121)
    ax.imshow(np.log(np.abs(tfrstft[:,:])))
    ax=fig.add_subplot(122)
    ax.imshow(np.log(np.abs(tfr[:,:])))
    plt.show()
