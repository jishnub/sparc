from __future__ import division
import numpy as np
import pyfits as pf
from matplotlib import pyplot as plt
import plotc
import scipy.ndimage.filters as filters
import scipy.integrate as scint
import fft_calculus as fc

Nsource=8
L=400*1e8
Rsun=695.9894*1e8




KBx=np.squeeze(pf.getdata("data/kernel/kernel_vectorpsi_01.fits"))
KBx=np.zeros_like(KBx)
KA=np.zeros_like(KBx)
KBz=np.zeros_like(KBx)
Kc=np.zeros_like(KBx)
Krho=np.zeros_like(KBx)
hess=np.zeros_like(KBx)

mag_field_dir="/home/jishnu/project/magnetic_inversions/working/start/"
Bz_start=np.squeeze(pf.getdata(mag_field_dir+"B_z2D.fits"))
Bx_start=np.squeeze(pf.getdata(mag_field_dir+"B_x2D.fits"))
rho_start=np.squeeze(pf.getdata(mag_field_dir+"density2D.fits"))
c_start=np.squeeze(pf.getdata(mag_field_dir+"soundspeed2D.fits"))

mag_field_dir="/home/jishnu/project/magnetic_inversions/working/true/"
Bz_true=np.squeeze(pf.getdata(mag_field_dir+"B_z2D.fits"))
Bx_true=np.squeeze(pf.getdata(mag_field_dir+"B_x2D.fits"))
rho_true=np.squeeze(pf.getdata(mag_field_dir+"density2D.fits"))
c_true=np.squeeze(pf.getdata(mag_field_dir+"soundspeed2D.fits"))


Nz,Nx=KBz.shape
area=Nx*Nz
x=np.linspace(0,L,Nx,endpoint=False)

polytropedir="/home/jishnu/project/magnetic_inversions/back/"
z=np.loadtxt(polytropedir+"polytrope")[:,0]*Rsun

x2d=np.atleast_2d(x)
z2d=np.atleast_2d(z).T
#~ A_trial=1e8*np.exp(-(x2d-L/2)**2/(L/20)**2)*np.exp(-(Rsun-z2d)/(Rsun/100))
#~ Bz_trial=fc.fft_deriv(A_trial,L=L)
dz=np.atleast_2d(np.gradient(z)).T
#~ Bx_trial=-np.gradient(A_trial)[0]/dz
#~ gl=plotc.gridlist(1,3)
#~ plotc.colorplot(A,subplot_index=next(gl))
#~ plotc.colorplot(Bz_trial,subplot_index=next(gl))
#~ plotc.colorplot(Bx_trial,subplot_index=next(gl))
#~ plt.show()
#~ exit()

kerneldir="/home/jishnu/project/magnetic_inversions/working/data/kernel/"

for n in xrange(1,Nsource+1):
	hess=np.squeeze(pf.getdata("data/kernel/hessian_0"+str(n)+".fits"))
	hess *= rho_start
	hess /= abs(hess).max()
 	hess[hess < 5e-3] = 5e-3
 	hess /= abs(hess).max() 	
	#~ hessn=np.squeeze(pf.getdata("data/kernel/hessian_0"+str(n)+".fits"))
	#~ hessn *= rho_start
	#~ hessn /= abs(hessn).max()
 	#~ hessn[hessn < 5e-3] = 5e-3
 	#~ hessn /= abs(hessn).max() 	
 	#~ hess+= hessn

	KA+=np.squeeze(pf.getdata(kerneldir+"kernel_vectorpsi_0"+str(n)+".fits"))/hess
	KBx+=np.squeeze(pf.getdata(kerneldir+"kernel_bx_0"+str(n)+".fits"))/hess
	KBz+=np.squeeze(pf.getdata(kerneldir+"kernel_bz_0"+str(n)+".fits"))/hess
	Kc+=np.squeeze(pf.getdata(kerneldir+"kernel_c_0"+str(n)+".fits"))/hess
	Krho+=np.squeeze(pf.getdata(kerneldir+"kernel_rho_0"+str(n)+".fits"))/hess

	#~ KA+=np.squeeze(pf.getdata(kerneldir+"kernel_vectorpsi_0"+str(n)+".fits"))
	#~ KBx+=np.squeeze(pf.getdata(kerneldir+"kernel_bx_0"+str(n)+".fits"))
	#~ KBz+=np.squeeze(pf.getdata(kerneldir+"kernel_bz_0"+str(n)+".fits"))
	#~ Kc+=np.squeeze(pf.getdata(kerneldir+"kernel_c_0"+str(n)+".fits"))
	#~ Krho+=np.squeeze(pf.getdata(kerneldir+"kernel_rho_0"+str(n)+".fits"))
	
#~ KA/=hess
#~ KBz/=hess
#~ KBx/=hess
#~ Kc/=hess
#~ Krho/=hess

f=np.fft.fftfreq(Nx)*Nx*2*np.pi/L*1j
f[Nx/2]=0
KBz_f=np.fft.fft(KBz)
f=np.atleast_2d(f)
dx_KBz=np.fft.ifft(f*KBz_f)
dx_KBz=np.real(dx_KBz)



dz_KBx=np.gradient(KBx)[0]/np.gradient(z)[0]

int_Bz_start=np.empty_like(Bz_start)
int_Bz_true=np.empty_like(Bz_true)
#~ int_Bz_trial=np.empty_like(Bz_trial)

for zind in xrange(Nz):
	int_Bz_start[zind]=fc.fft_int(Bz_start[zind],L)
	int_Bz_true[zind]=fc.fft_int(Bz_true[zind],L)
	#~ int_Bz_trial[zind]=fc.fft_int(Bz_trial[zind],L)

int_Bx_start=scint.cumtrapz(Bx_start[:,0],x=z,initial=0)
int_Bx_true=scint.cumtrapz(Bx_true[:,0],x=z,initial=0)
#~ int_Bx_trial=scint.cumtrapz(Bx_trial[:,0],x=z,initial=0)

#~ int_Bx_start,int_Bx_true,int_Bx_trial=np.atleast_2d(int_Bx_start,int_Bx_true,int_Bx_trial)
int_Bx_start,int_Bx_true=np.atleast_2d(int_Bx_start,int_Bx_true)
int_Bx_start=int_Bx_start.T
int_Bx_true=int_Bx_true.T
#~ int_Bx_trial=int_Bx_trial.T

A_start=int_Bz_start - int_Bx_start
A_true=int_Bz_true - int_Bx_true
#~ A_trial=int_Bz_trial - int_Bx_trial
A_start=A_start/Rsun
A_true=A_true/Rsun
#~ A_trial=A_trial/Rsun

z=z-Rsun
yax_ll=z[200]
yax_ul=None

#~ KBz=filters.gaussian_filter1d(KBz,  2)
#~ KBx=filters.gaussian_filter1d(KBx,  2)
#~ Kc=filters.gaussian_filter1d(Kc,  2)
#~ Krho=filters.gaussian_filter1d(Krho,2)
#~ KA=filters.gaussian_filter1d(KA,    2)
#~ Kw=filters.gaussian_filter1d(Kw,    2)

Kc=filters.gaussian_filter(Kc,  [2,2])
Krho=filters.gaussian_filter(Krho,[4,2])
KA=filters.gaussian_filter(KA,    [2,2])
KBx=filters.gaussian_filter(KBx,  [2,2])
KBz=filters.gaussian_filter(KBz,  [2,2])

inner_cutoff = 220
surf_cutoff=-1
nrows=3;ncols=2
gridlist=plotc.gridlist(nrows,ncols)

#~ plotc.colorplot(KBx,y=z,subplot_index=next(gridlist),yr=[z[inner_cutoff],yax_ul],
#~ title="KBx",ylabel="Depth from surface (Mm)")
#~ plotc.colorplot(KBz,y=z,subplot_index=next(gridlist),yr=[z[inner_cutoff],yax_ul],
#~ title="KBz",ylabel="Depth from surface (Mm)")

KA_inner_cutoff=270

KA_plot=KA[KA_inner_cutoff:]
plotc.colorplot(KA_plot,y=z[KA_inner_cutoff:],subplot_index=next(gridlist),title="A Kernel",vmin=-KA_plot.max())
dA=A_true-A_start;dA=dA[KA_inner_cutoff:]
plotc.colorplot(dA,y=z[KA_inner_cutoff:],subplot_index=next(gridlist),title="A true - A start")

Kcplot=Kc[inner_cutoff:surf_cutoff]
plotc.colorplot(Kcplot,y=z[inner_cutoff:surf_cutoff],subplot_index=next(gridlist),cbar_sci=True,
title="c Kernel",vmin=-Kcplot.max())
dc=(c_true[inner_cutoff:surf_cutoff]-c_start[inner_cutoff:surf_cutoff])/c_start[inner_cutoff:surf_cutoff]
plotc.colorplot(dc,y=z[inner_cutoff:surf_cutoff],subplot_index=next(gridlist),cbar_sci=True,
title="c true - c start")


Krhoplot=Krho[inner_cutoff:surf_cutoff]
plotc.colorplot(Krhoplot,y=z[inner_cutoff:surf_cutoff],subplot_index=next(gridlist),cbar_sci=True,
title="rho Kernel",vmin=-Krhoplot.max())
drho=(rho_true[inner_cutoff:surf_cutoff]-rho_start[inner_cutoff:surf_cutoff])/rho_start[inner_cutoff:surf_cutoff]
plotc.colorplot(drho,y=z[inner_cutoff:surf_cutoff],subplot_index=next(gridlist),cbar_sci=True,
title="rho true - rho start")

plt.tight_layout()
plt.show()

write_vector_potential=False
if write_vector_potential: 
	pf.writeto('start/model_vectorpsi.fits',A_start,clobber=True)
	pf.writeto('data/model_vectorpsi.fits',A_start,clobber=True)
	pf.writeto('data/update/model_vectorpsi_00.fits',A_start,clobber=True)
	#~ pf.writeto('model_vectorpsi.fits',A_start,clobber=True)
	pf.writeto('field/true_vectorpsi.fits',A_true,clobber=True)
	pf.writeto('true/true_vectorpsi.fits',A_true,clobber=True)



