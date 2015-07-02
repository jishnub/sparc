import numpy as np
import plotc
from matplotlib import pyplot as plt
import pyfits as pf

Btruedir="/home/jishnu/project/magnetic_inversions/working/true/"

z=np.loadtxt('/home/jishnu/project/magnetic_inversions/back/polytrope')[:,0];Rsun=695.8e3;z=(z-z[-1])*Rsun
x=np.linspace(-200,200,256,endpoint=False)

Bztrue=pf.getdata(Btruedir+'B_z2D.fits')
Bxtrue=pf.getdata(Btruedir+'B_x2D.fits')

soundspeed=pf.getdata(Btruedir+'soundspeed2D.fits')
density=pf.getdata(Btruedir+'density2D.fits')

calf=np.sqrt((Bztrue**2+Bxtrue**2)/density/(4*np.pi))
cratio=calf/soundspeed

plotc.quiver2D(Bxtrue,Bztrue,x=x,y=z,every=[3,1],sp=121,yr=[-0.5e3,None],xr=[-100,100],key=True,
keysuffix=" Gauss",title="Magnetic Field",titley=1.01,usetex=True,
xlabel="Horizontal Distance (Mm)",ylabel="Depth (km)")
plotc.colorplot(cratio,x=x,y=z,sp=122,yr=[-0.5e3,None],xr=[-100,100],vmax=2,title="$c_A/c_S$",titley=1.01,usetex=True,
xlabel="Horizontal Distance (Mm)")

plt.setp(plt.gca().get_yticklabels(),visible=False)
plt.subplots_adjust(wspace=0)

plotc.figuresize(13,5)
plotc.show()
#~ plotc.savefig('/home/jishnu/project/departmental_talk/field_cratio.eps')
