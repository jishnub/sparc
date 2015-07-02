from matplotlib import pyplot as plt,cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


Xax,Yax=np.meshgrid(np.linspace(-1,1,15),np.linspace(-1,1,15))
z = Xax**2+10*Xax**2*Yax**2+Yax**2-Xax**4*Yax**4

ax=plt.subplot(121, projection='3d')
ax.plot_surface(Xax, Yax, z,  rstride=1, cstride=1,shade=False,cmap='coolwarm',linewidth=0.1)
ax.scatter(Xax,Yax,z,color='Brown')
plt.xlim(-1,1)
plt.ylim(-1,1)
t1=plt.title("Forward problem")
t1.set_y(0.95)

ax=plt.subplot(122, projection='3d')
ax.plot_surface(Xax, Yax, z,  rstride=1, cstride=1,cmap='coolwarm',linewidth=0.1)
x=np.linspace(-0.9,-0.1,6);y=np.linspace(0.9,0.1,6);zp=x**2+10*x**2*y**2+y**2-x**4*y**4
ax.plot(x,y,zp,'o',color='orange',markersize=7)
ax.plot(x,y,zp,'-',color='Brown')
t2=plt.title("Inverse Problem")
t2.set_y(0.95)

plt.gcf().set_size_inches(15,5)
#~ plt.show()
plt.savefig('../presentation/forawrd_inverse.png')


