import numpy as np
import plotc
from matplotlib import pyplot as plt
import pyfits as pf

truedir="/home/jishnu/project/magnetic_inversions/working/true/"
polytropefile='/home/jishnu/project/magnetic_inversions/back/polytrope_sh'

z,_,densitypoly,pressurepoly,_,_=np.loadtxt(polytropefile,unpack=True);Rsun=695.8e3;z=(z-z[-1])*Rsun
x=np.linspace(-200,200,256,endpoint=False)

density=pf.getdata(truedir+'density2D.fits')
pressure=pf.getdata(truedir+'pressure2D.fits')

densitypoly2D=np.atleast_2d(densitypoly).T
densitypoly2D=np.tile(densitypoly2D,[1,256])

drho=density-densitypoly2D

pressurepoly2D=np.atleast_2d(pressurepoly).T
pressurepoly2D=np.tile(pressurepoly2D,[1,256])

dp=pressure-pressurepoly2D

plotc.colorplot(dp/pressurepoly2D,x=x,y=z,sp=121,yr=[-0.2e3,None],xr=[-100,100],
xlabel="Horizontal Distance (Mm)",ylabel="Depth (km)",usetex=True,
title="Relative Pressure Difference",titley=1.02)
plotc.colorplot(drho/densitypoly2D,x=x,y=z,sp=122,yr=[-0.2e3,None],xr=[-100,100],
xlabel="Horizontal Distance (Mm)",ylabel="Depth (km)",usetex=True,
title="Relative Density Difference",titley=1.02)

#~ plotc.show()
plotc.figuresize(13,5)
plotc.savefig('/home/jishnu/project/departmental_talk/dp_drho.eps')
