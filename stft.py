import numpy as np

def tfrstft(x,t=None,N=None,f=None,itc=None,h=None,nwindow=4,silent=0,use_nufft=False):
    xrow = np.shape(x)[0] 
    if t is None: t=np.array(range(0,xrow),dtype=np.int)
    if N is None: N=xrow 

    if itc is None:  
        Nt=N
        itc=range(0,Nt)
    else:
        Nt=len(itc)

    if h is None:
        hlength=np.int(np.floor(N/nwindow))
        hlength=hlength+1-np.mod(hlength,2)
        h=0.54 - 0.46*np.cos(2.0*np.pi*np.array(range(0,hlength))/(hlength+1)) #Hamming
    h=h/np.linalg.norm(h)
    hrow=len(h)
    Lh=np.int(np.round((hrow-1)/2)) ##??

    tfr=np.zeros((N,Nt),dtype=np.complex) # plane by default
    ti=t[itc]
    Np2=np.int(N/2)
    for icol in range(0,Nt):       
        ti=t[itc[icol]]
        taumin=np.int(-np.min([np.floor(Np2)-1,Lh,ti]))
        taumax=np.int(np.min([np.floor(Np2)-1,Lh,xrow-ti-1]))
        i0=np.int(np.mod(N+taumin+1,N))
        i1=np.int(np.mod(N+taumax+1,N))
        tfr[ti+taumin:ti+taumax,icol]=x[ti+taumin:ti+taumax]*(h[Lh+taumin:Lh+taumax])
    ###Choose FFT or DFT
    if f is None:
        for i in range(0,Nt):
            tfr[:,i]=np.fft.fft(tfr[:,i])
        return tfr
    elif use_nufft:
        print("Nufft is not implemented yet.")
        return tfr
if __name__ == "__main__":
    import sampledata as sd
    import matplotlib.pyplot as plt
    import time
    start=time.time()
    nsamp=1024*4
    t,x=sd.genmultifm623(nsamp)
    stftx=tfrstft(x)
    print(time.time()-start,"sec")
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.imshow(np.abs(stftx[:,:]))
    plt.show()
    
