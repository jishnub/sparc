import numpy as np
import plotc
import pyfits as pf

it=0

graddir="/home/jishnu/project/magnetic_inversions/working/data/update/"
grad=pf.getdata(graddir+"gradient_vectorpsi_0"+str(it)+".fits")
grad/=grad.max()
z=np.loadtxt('/home/jishnu/project/magnetic_inversions/back/polytrope')[:,0];Rsun=695.8e3;z=(z-z[-1])*Rsun
x=np.linspace(-200,200,256,endpoint=False)

Atruedir="/home/jishnu/project/magnetic_inversions/working/true/"
Astartdir="/home/jishnu/project/magnetic_inversions/working/start/"

Atrue=pf.getdata(Atruedir+'true_vectorpsi.fits')
Astart=pf.getdata(Astartdir+'model_vectorpsi.fits')
Adiff=Atrue-Astart

plotc.colorplot(Adiff,x=x,y=z,yr=[-0.8e3,None],xr=[-100,100],vmax=1,centerzero=True,xbins=5,
title=r"$\Delta A$",xlabel="Horizontal Distance (Mm)",ylabel="Depth (km)",usetex=True,y_sci=False
,sp=121)
plotc.colorplot(grad,x=x,y=z,yr=[-0.8e3,None],xr=[-100,100],vmax=1,centerzero=True,xbins=5,
title=r"Normalized Vector Potential Kernel",xlabel="Horizontal Distance (Mm)",ylabel="Depth (km)",usetex=True,y_sci=False
,sp=122)
plotc.figuresize(13,5)
plotc.tight_layout()
plotc.savefig('/home/jishnu/project/magnetic_inversions/report/vector_potential_kernel_'+str(it)+'.eps')
#~ plotc.show()
