from __future__ import division
import numpy as np
import sys,os
import plotc
import pyfits
import deriv6ord,fftcalculus as fc
import warnings

Rsun=695.9894
zeroB=True

def derivx(arrxz): return fc.differentiate(arrxz,axis=0,period=Lx)

def derivz(arrxz):
    darr=np.asfortranarray(np.zeros(arrxz.shape,dtype=float))
    deriv6ord.dbyd2(darr,arrxz,ibc=1)
    return darr/dz

codedir=os.path.join(os.environ['HOME'],"sparc")
modelp = np.loadtxt(os.path.join(codedir,'polytrope'))
user=os.environ['USER']
modeldir=os.path.join('/scratch',user,'magnetic')

nx = 256
Lx = 400
cutoff = np.floor(nx/3)
x = np.linspace(-Lx/2,Lx/2,num=nx,endpoint=False)


z = (modelp[:,0]-1.)*Rsun;
nz=z.shape[0]
z2d=np.asfortranarray(np.atleast_2d(z))
dz2d=np.asfortranarray(np.zeros_like(z2d))
deriv6ord.dbyd2(dz2d,z2d,ibc=1)
dz=np.squeeze(dz2d)

p0 = modelp[:,3]
rho0 = modelp[:,2]
c2 = modelp[:,1]**2

grav=2.775e4
convert = 1e8



A0=5000
decay_x=10;decay_z=10;

fx=-(x/decay_x)/(1+(x/decay_x)**4)
dfx=-decay_x**3*(decay_x**4-3*x**4)/(decay_x**4+x**4)**2
gz=np.exp(z/decay_z)
dgz=1/decay_z*gz


bx=-A0*np.atleast_2d(fx).T*dgz
bz=A0*np.atleast_2d(dfx).T*gz

A_true=A0*np.atleast_2d(fx).T*gz

decay_z_start=12
gz_start = np.exp(z/decay_z_start)
A_start=A0*np.atleast_2d(fx).T*gz_start

bx_start=-derivz(A_start)
bz_start=fc.differentiate(A_start,axis=0,period=Lx)

if zeroB:
    bx=np.zeros_like(bx)
    bz=np.zeros_like(bz)
    bx_start=np.zeros_like(bx)
    bz_start=np.zeros_like(bx)
    A_start=np.zeros_like(A_start)
    A_true=np.zeros_like(A_true)

force = derivz(bx)-derivx(bz);

hz = bx*force;
hx = bz*force;

p=fc.integrate(hx,axis=0,period=Lx)

p-= p[0] - p0


dzp = derivz(p);

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
soundspeed=np.sqrt(soundspeed_squared)
alfven_speed=np.sqrt((bx**2+bz**2)/rho)
alf_c_ratio=alfven_speed/soundspeed

plotc.plt.figure()
plotc.colorplot(bz.T*np.sqrt(4*np.pi),y=z,x=x,sp=221,title="Bz true")
plotc.colorplot(bz_start.T*np.sqrt(4*np.pi),y=z,x=x,sp=222,title="Bz start")
plotc.colorplot(bx.T*np.sqrt(4*np.pi),y=z,x=x,sp=223,title="Bx true")
plotc.colorplot(bx_start.T*np.sqrt(4*np.pi),y=z,x=x,sp=224,title="Bx start")


#~ plotc.plt.figure()
#~ plotc.colorplot(dp.T,x=x,y=z,title="Relative Pressure Difference",sp=121,
                #~ ylim=[-2.5,-0.39],
                #~ axes_properties=dict(xlabel="x",ylabel="z",
                #~ locator_properties_x=dict(nbins=4),
                #~ locator_properties_y=dict(nbins=5)),
                #~ colorbar_properties=dict(orientation='horizontal'
                #~ ,shrink=0.9))
#~ 
#~ plotc.colorplot(drho.T,x=x,y=z,title="Relative Density Difference",sp=122,
                #~ ylim=[-2.5,-0.39],
                #~ axes_properties=dict(xlabel="x",
                #~ locator_properties_x=dict(nbins=4),
                #~ locator_properties_y=dict(nbins=5),
                #~ hide_yticklabels=True),
                #~ centerzero=True,
                #~ colorbar_properties=dict(orientation='horizontal',shrink=0.9))
#~ 
#~ plotc.plt.subplots_adjust(wspace=0)

#~ plotc.plt.figure()

#~ plotc.quiver2D(bx.T,bz.T,x=x,y=z,every=[2,15],xr=[-50,50],scale=5000,key=True)
#~ plotc.plt.show()

write_models=True

if write_models:
    if not os.path.exists(os.path.join(modeldir,'true')): os.makedirs(os.path.join(modeldir,'true'))
    if not os.path.exists(os.path.join(modeldir,'start')): os.makedirs(os.path.join(modeldir,'start'))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pyfits.writeto(os.path.join(modeldir,'true','B_z2D.fits'),bz.T*np.sqrt(4*np.pi),clobber=True)
        pyfits.writeto(os.path.join(modeldir,'true','B_x2D.fits'),bx.T*np.sqrt(4*np.pi),clobber=True)
        pyfits.writeto(os.path.join(modeldir,'true','true_vectorpsi.fits'),A_true.T*np.sqrt(4*np.pi)/Rsun,clobber=True)
        pyfits.writeto(os.path.join(modeldir,'start','model_vectorpsi_ls00.fits'),A_start.T*np.sqrt(4*np.pi)/Rsun,clobber=True)
        
        pyfits.writeto(os.path.join(modeldir,'true','pressure2D.fits'),p.T,clobber=True)
        pyfits.writeto(os.path.join(modeldir,'start','pressure2D.fits'),p.T,clobber=True)
        pyfits.writeto(os.path.join(modeldir,'true','density2D.fits'),rho.T,clobber=True)
        pyfits.writeto(os.path.join(modeldir,'start','density2D.fits'),rho.T,clobber=True)
        pyfits.writeto(os.path.join(modeldir,'true','soundspeed2D.fits'),soundspeed.T,clobber=True)
        pyfits.writeto(os.path.join(modeldir,'start','soundspeed2D.fits'),soundspeed.T,clobber=True)
        pyfits.writeto(os.path.join(modeldir,'start','model_c_ls00.fits'),soundspeed.T,clobber=True)


