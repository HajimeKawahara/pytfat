if __name__ == "__main__":

    import smethod 
    import stft
    import sys
    import numpy as np
    import matplotlib.pyplot as plt
    import time
    start=time.time()
    
    dat=np.loadtxt("KRD_SFHx_cross.dat")
    x=dat[2950:2950+1024*4]
    print(len(x))
    

    tfr=smethod.tfrsm(x,Lp=12,nwindow=4)
#    tfrstft=stft.tfrstft(x,nwindow=8)
    print(time.time()-start,"sec")
    
    fig=plt.figure()
    ax=fig.add_subplot(121)
#    c=ax.imshow(np.log10((np.abs(tfrstft[0:256,:])**2)),cmap="CMRmap",vmin=-44,vmax=-41.5)
#    plt.colorbar(c)
    ax.set_aspect(0.7/ax.get_data_ratio())
    plt.gca().invert_yaxis()

    
    ax=fig.add_subplot(122)
    c=ax.imshow(np.log10((np.abs(tfr[0:256,:]))),cmap="CMRmap",vmin=-44,vmax=-41.5)
    plt.colorbar(c)
    ax.set_aspect(0.7/ax.get_data_ratio())
    plt.gca().invert_yaxis()
    plt.show()
