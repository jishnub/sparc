from __future__ import division
import numpy as np
import math
import sys
import plotc
import pyfits

def derivx(eyekx2D,arrxz):
    darrxz=np.fft.fft(bz,axis=0)
    darrxz*=eyekx2D
    darrxz=np.real(np.fft.ifft(darrxz,axis=0))
    return darrxz

def derivz(z,arrxz):
    darrxz=np.gradient(arrxz)[1]
    dz=np.gradient(z)
    return darrxz/dz

modelp = np.loadtxt('polytrope')
nx = 256
Lx = 400
cutoff = math.floor(nx/3)
x = np.linspace(-Lx/2,Lx/2,num=nx,endpoint=False)
dx=x[1]-x[0]
x2D=np.atleast_2d(x).T
kx = np.fft.fftfreq(nx,dx)*2*np.pi
eyekx = 1j*kx;
eyekx2D = np.atleast_2d(eyekx).T
z = (modelp[:,0]-1.)*695.9894;
zfull=z;
zmin = -5.2;
zmax = z.max();
#~ nz1 = np.argmin(abs(z-zmin))
nz1=0
nz2 = len(z)
nz = nz2-nz1
z = z[nz1:nz2]
p0 = modelp[nz1:nz2,3]
rho0 = modelp[nz1:nz2,2]
c2 = modelp[nz1:nz2,1]**2
print c2.shape
#~ grav = modelp[nz1:nz2,3];

grav=2.775e4
convert = 1e8

zcut = nz-1
sigma = 8. + 14./(1. + np.exp((z[zcut]-z)/5))
sigma2D=np.atleast_2d(sigma)

psi = np.zeros((nx,nz))
temp = np.exp(np.polyval(np.polyfit(z,np.log(p0),3),z));
temp2D = np.atleast_2d(temp)
p0fn= np.polyfit(z,p0,4)
zpts= np.arange(-5,0)
p0pts=np.polyval(p0fn,zpts)
p0pts[1]=6.5e+7
p0pts[0]=6e+7
coeff=np.polyfit(zpts,p0pts,4)

bz=temp2D**0.5*(1/(1. + np.exp((x2D-sigma2D)/2.5)) -1/(1.+ np.exp((x2D+sigma2D)/2.5)))
fluxes=np.atleast_2d(np.sum(bz,axis=0)*dx)
bz/=fluxes

normalization = 500/bz[:,zcut].max()
bz*=normalization

psi=-np.fft.fft(bz,axis=0)
psi[0]=0 
psi[1:]/=eyekx2D[1:]
psi=np.real(np.fft.ifft(psi,axis=0))

dzbz=derivz(z,bz);
dzbz=dzbz[:,zcut];
bx = derivz(z,psi);

force = -derivx(eyekx2D,bz)+ derivz(z,bx);

hz = bx*force;
hx = bz*force;

p=np.fft.fft(hx,axis=0)
p[0]=0 
p[1:]/=eyekx2D[1:]
p=np.real(np.fft.ifft(p,axis=0))
p-= p[0] - p0


dzp = derivz(z,p);

dp = p/p0 - 1

if dp.min() < -1 :
    print "Negative pressure"
    sys.exit(0)


rho = -(hz + dzp)/(grav*convert);
rho-=rho[0]- rho0
drho=rho/rho0-1



if drho.min() < -1 :
    print "Negative density"
    sys.exit(0)
    


soundspeed_squared=c2*rho0/p0*p/rho # This works even if c2 etc are 1D, since x is the leading dimension and multiplication is along z
soundspeed=np.sqrt(c2mod)
alfven_speed=np.sqrt((bx**2+bz**2)/rho)
alf_c_ratio=alfven_speed/soundspeed

plotc.colorplot(dp.T,x=x,y=z,title="Relative Pressure Difference",sp=121,
                ylim=[-2.5,-0.39],
                axes_properties=dict(xlabel="x",ylabel="z",
                locator_properties_x=dict(nbins=4),
                locator_properties_y=dict(nbins=5)),
                colorbar_properties=dict(orientation='horizontal',
                ticks=[-0.25*i for i in xrange(1,4)],shrink=0.9))
#~ 
plotc.colorplot(drho.T,x=x,y=z,title="Relative Density Difference",sp=122,
                ylim=[-2.5,-0.39],
                axes_properties=dict(xlabel="x",
                locator_properties_x=dict(nbins=4),
                locator_properties_y=dict(nbins=5),
                hide_yticklabels=True),
                centerzero=True,
                colorbar_properties=dict(orientation='horizontal',shrink=0.9))

plotc.plt.subplots_adjust(wspace=0)

plotc.plt.figure()

plotc.quiver2D(bx.T,bz.T,x=x,y=z,every=[2,15],xr=[-50,50])
plotc.plt.show()

pyfits.writeto('true/B_z2D.fits',bx.T*np.sqrt(4*np.pi))
pyfits.writeto('true/B_x2D.fits',bx.T*np.sqrt(4*np.pi))
pyfits.writeto('true/pressure2D.fits',p.T)
pyfits.writeto('true/density2D',rho.T)
pyfits.writeto('true/soundspeed2D',soundspeed.T)


