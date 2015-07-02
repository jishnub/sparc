from __future__ import division
import numpy as np
import pyfits as pf
import plotc
import scipy.integrate
import fftcalculus as fc
import os


Nsource=8
L=400*1e8
Rsun=695.9894*1e8
L=L/Rsun

codedir="/home/jishnu/sparc"
start_field_dir="/scratch/jishnu/magnetic/start/"
#~ Bz_start=np.squeeze(pf.getdata(os.path.join(start_field_dir,"B_z2D.fits")))
#~ Bx_start=np.squeeze(pf.getdata(os.path.join(start_field_dir,"B_x2D.fits")))

true_field_dir="/scratch/jishnu/magnetic/true/"
Bz_true=np.squeeze(pf.getdata(os.path.join(true_field_dir,"B_z2D.fits")))
Bx_true=np.squeeze(pf.getdata(os.path.join(true_field_dir,"B_x2D.fits")))

Nz,Nx=Bz_true.shape
area=Nx*Nz
x=np.linspace(0,L,Nx,endpoint=False)
z=np.loadtxt(os.path.join(codedir,"polytrope"))[:,0]


Bz_start=Bz_true*

int_Bz_start=np.empty_like(Bz_start)
int_Bz_true=np.empty_like(Bz_true)

int_Bz_start=fc.integrate(Bz_start,axis=1,period=L)
int_Bz_true=fc.integrate(Bz_true,axis=1,period=L)

int_Bx_start=scipy.integrate.cumtrapz(Bx_start[:,0],x=z,initial=0)
int_Bx_true=scipy.integrate.cumtrapz(Bx_true[:,0],x=z,initial=0)

int_Bx_start,int_Bx_true=np.atleast_2d(int_Bx_start,int_Bx_true)
int_Bx_start=int_Bx_start.T
int_Bx_true=int_Bx_true.T

A_start=int_Bz_start - int_Bx_start
A_true=int_Bz_true - int_Bx_true

z-=1
yax_ll=z[200]
yax_ul=None
inner_cutoff=270
surf_cutoff=-1
gridlist=plotc.gridlist(1,2)

ploty=z[inner_cutoff:]*Rsun/1e8

dA=A_true-A_start;dA=dA[inner_cutoff:]
plotc.colorplot(A_start[inner_cutoff:],y=ploty,subplot_index=next(gridlist),title="A start")
plotc.colorplot(dA,y=ploty,subplot_index=next(gridlist),title="A true - A start")

plotc.plt.show()
write_vector_potential=False
if write_vector_potential: 
	pf.writeto(os.path.join(start_field_dir,'model_vectorpsi_ls00.fits'),A_start,clobber=True)
	pf.writeto(os.path.join(true_field_dir,'true_vectorpsi.fits'),A_true,clobber=True)
