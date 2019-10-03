import numpy as np
import scipy 
import stft 



def compute_crossconj_all(x,y,N,Lp):

    arr=[]
    for k in range(0,N):
        i=np.max([k-Lp,0])
        j=np.min([k+Lp+1,N])
        arr.append(np.sum(np.real(x[:,i:j]*np.conj(y[:,i:j][:,::-1])),axis=1))
    arr=np.array(arr).T
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

    sm=compute_crossconj_all(tfrstftx,tfrstfty,nsamplef,Lp)

    return sm

if __name__ == "__main__":
    import sampledata as sd
    import matplotlib.pyplot as plt
    import time
    start=time.time()
    nsamp=512
    t,x=sd.genmultifm622(nsamp)
    tfr=tfrsm(x,Lp=32,nwindow=8)
    tfrstft=stft.tfrstft(x,nwindow=8)
    print(time.time()-start,"sec")
    
    fig=plt.figure()
    ax=fig.add_subplot(121)
    ax.imshow((np.abs(tfrstft[:,:])))
    ax=fig.add_subplot(122)
    ax.imshow((np.abs(tfr[:,:])))
    plt.show()
    
