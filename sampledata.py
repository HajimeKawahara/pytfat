import numpy as np

def genfm(nsamp,wp=1.0,wf=0.01,s=1.0,xend=100*np.pi):
    x=np.linspace(0.0,xend,nsamp)
    y=np.cos(wp*x+s*np.sin(wf*x))
    iw=wp+s*wf*np.cos(wf*x)
    ynorm=np.pi/x[-1]
    return x, y, iw, ynorm

def genfmX(nsamp,Fp=1.0,Ff=0.01,amp=0.01,xend=100*np.pi):
    #amp [fm fraction]
    x=linspace(0.0,xend,nsamp)
    y=np.cos(2.0*np.pi*Fp*x+amp*Fp/Ff*np.sin(2.0*np.pi*Ff*x))
    iF=Fp+amp*Fp*np.cos(2.0*np.pi*Ff*x)
    ynorm=pi/x[-1]
    return x, y, iF, ynorm

def genlinfm(nsamp,wp=1.0,s=0.01, x=100*np.pi):
    x=np.linspace(0.0,x,nsamp)
    y=np.cos(wp*x+0.5*s*x*x)
    iw=wp+s*x
    ynorm=np.pi/x[-1]
    return x, y, iw, ynorm


def genstepfm(nsamp=1024,x=1.0):
    x=np.linspace(0.0,x,nsamp)
    iw=nsamp/4*np.atan(250*(x-0.5))+np.pi*nsamp/4
    phase=cumsum(iw)/nsamp
    y=np.cos(phase)
    iy=np.sin(phase)*im
    ynorm=np.pi/x[-1]
    return x, y+iy, iw, ynorm


def genmultifm622(nsamp=512):
    #Boashash+2015 Example 6.2.2
    t=(np.linspace(-1.0,1.0,nsamp))
    x=np.exp(-t*t)*(np.cos(25.0*np.pi*t))+np.cos(120.0*t*t*t+45.0*np.pi*t)+1.5*np.exp(-25.0*t*t)*np.cos(40.0*np.pi*t*t+150.0*np.pi*t)
    return t, x


def genmultifm622x(nsamp=512):
    #Boashash+2015 modified Example 6.2.2
    t=(np.linspace(-1.0,1.0,nsamp))
    x=np.exp(-t*t)*np.cos(25.0*np.pi*t)+np.cos(12.0*t*t*t+40.0*np.pi*t)+1.5*np.exp(-25.0*t*t)*np.cos(4.0*np.pi*t*t+65.0*np.pi*t)
    return t, x


def genmultifm623(nsamp=256):
    #Boashash+2015 Example 6.2.3
    t=(np.linspace(-1.0,1.0,nsamp))
    x=np.cos(20.0*np.sin(np.pi*t)+30.0*np.pi*t)+np.sin(20.0*np.cos(np.pi*t)+100.0*np.pi*t)
    return t, x


def genmultifm623m(nsamp=256):
    t=(np.linspace(-1.0,1.0,nsamp))
    x=np.cos(40.0*np.pi*t)+np.sin(20.0*np.cos(np.pi*t)+100.0*np.pi*t)
    return t, x



def genmultifm623x(nsamp=256):
    #Boashash+2015 Example 6.2.3 modified
    t=(np.linspace(-1.0,1.0,nsamp))
    x=np.cos(2.0*np.sin(np.pi*t)+30.0*np.pi*t)+np.sin(2.0*np.cos(np.pi*t)+60.0*np.pi*t)
    return t, x

