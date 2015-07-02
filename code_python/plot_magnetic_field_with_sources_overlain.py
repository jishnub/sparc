import numpy as np
import plotc
from matplotlib import pyplot as plt
import pyfits as pf

Btruedir="/home/jishnu/project/magnetic_inversions/working/true/"
z=np.loadtxt('/home/jishnu/project/magnetic_inversions/back/polytrope')[:,0];Rsun=695.8e3;z=(z-z[-1])*Rsun
x=np.linspace(-200,200,256,endpoint=False)

Bztrue=pf.getdata(Btruedir+'B_z2D.fits')
Bxtrue=pf.getdata(Btruedir+'B_x2D.fits')

plotc.quiver2D(Bxtrue,Bztrue,x=x,y=z,every=[3,2],yr=[-0.2e3,None],xr=[-100,100],key=True,
keysuffix=" Gauss",title="Magnetic Field",titley=1.01)

masterpix='/home/jishnu/project/magnetic_inversions/working/data/master.pixels'
sourcelocs=np.loadtxt(masterpix)

circ1=plt.Circle((sourcelocs[0],-1),radius=3,color='#AA2222')
circ2=plt.Circle((sourcelocs[7],-1),radius=3,color='#AA2222')
plt.gca().add_patch(circ1)
plt.gca().add_patch(circ2)

R1=10

plotc.show()
