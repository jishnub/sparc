from __future__ import division
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import patches
import pyfits as pf
import plotc

Rsun=695.8
dt = 30

pi=np.pi
cos=np.cos
Poly=np.zeros(3)
Poly[0]=0.0010
Poly[1]=0.0030
Poly[2]=-0.0006
f_mode_const=0.00293795
df = 0.5
f_low = 1.1

source=8
master_pixels=np.loadtxt('/home/jishnu/project/magnetic_inversions/working/data/master.pixels')
source=str(source).zfill(2)
vzcc=pf.getdata('/home/jishnu/project/magnetic_inversions/working/quiet.fits')
#~ vzcc=pf.getdata('/home/jishnu/project/magnetic_inversions/working/data/forward'+source+'/vz_cc.fits')
#~ vzcc=pf.getdata('/home/jishnu/project/magnetic_inversions/working/vz_unfiltered.fits')
#~ vzcc=pf.getdata('/home/jishnu/project/magnetic_inversions/working/fmode.fits')
#~ dat=pf.getdata('/home/jishnu/project/magnetic_inversions/working/data/forward'+source+'/data.fits')
vzcc=np.squeeze(vzcc)
nxt=vzcc.shape
vzkw=np.fft.fft2(vzcc)
#~ datkw=np.fft.fft2(dat)
powsp=abs(vzkw)**2

nx=nxt[1]
nt=nxt[0]

fmode=np.zeros((nx,1,nt))
xlength=400*10**8
L=xlength*10.**(-8.)*nx/(nx-1.)
dx=L/nx

kf = np.fft.fftfreq(nx,dx)*2*pi
ff = np.fft.fftfreq(nt,dt)
k = abs(kf)
w = abs(ff) * 2.*pi

f0=0.7*f_mode_const*k**0.5*1e3
f1=(Poly[0] + Poly[1]*k +Poly[2]*k**2.)*1e3
f = w/(2.*pi)*1e3
ff=ff*1e3

kax,fax=np.meshgrid(k,f)

x=np.linspace(-xlength/2,xlength/2,nx)/1e8
t=np.linspace(0,nt*dt/60,nt)
Xax,Tax=np.meshgrid(x,t)

for i in xrange(nx):
  delta = (f1[i] - f0[i])
  for j in xrange(nt):
    d = f[j] - f0[i]  
    if d>0 and d < delta:
      fmode[i,0,j] = 0.5*(1.+cos(pi*(d-delta*0.5)/(delta*0.5)))


for j in xrange(nt):
  if (f[j] < f_low+df):
    fmode[:,0,j] = fmode[:,0,j]*0.5*(1.+cos(pi*(f[j]-(f_low+df))/df) )
  if (f[j] < f_low): fmode[:,0,j] = 0.


fmode=np.squeeze(fmode).T


g=plotc.gridlist(1,2)

good_quadrant=powsp[:nt/2,:nx/2]
goodk=kf[:nx/2]*Rsun
goodf=ff[:nt/2]
plotc.colorplot(good_quadrant,x=goodk,y=goodf,subplot_index=next(g),
xlabel="k$R_\odot$ ",ylabel="$\omega$ (mHz)",title="Power spectrum",cbar_sci=True,
yr=[0,5],xr=[0,1.5*Rsun],usetex=True,x_sci=False)

plotc.gca().text(550,2,"f-mode",fontsize="15")
plotc.gca().text(550,3.4,"p1-mode",fontsize="15")

#~ filter_low_cutoff=np.where(f0>f_low)[0][0]
#~ k_filt_line_low=k[filter_low_cutoff-1:nx/2]*Rsun
#~ f0_filt_line_low=f0[filter_low_cutoff-1:nx/2]
#~ plt.plot(k_filt_line_low,f0_filt_line_low,zorder=2,color='brown')
#~ plt.plot(-k_filt_line_low,f0_filt_line_low,zorder=2,color='brown')
#~ plt.plot(k_filt_line_low,-f0_filt_line_low,zorder=2,color='brown')
#~ plt.plot(-k_filt_line_low,-f0_filt_line_low,zorder=2,color='brown')
#~ 
#~ filter_low_cutoff=np.where(f1>f_low)[0][0]
#~ k_filt_line_low=k[filter_low_cutoff-1:nx/2]*Rsun
#~ f1_filt_line_low=f1[filter_low_cutoff-1:nx/2]
#~ plt.plot(k_filt_line_low,f1_filt_line_low,zorder=2,color='brown')
#~ plt.plot(-k_filt_line_low,f1_filt_line_low,zorder=2,color='brown')
#~ plt.plot(k_filt_line_low,-f1_filt_line_low,zorder=2,color='brown')
#~ plt.plot(-k_filt_line_low,-f1_filt_line_low,zorder=2,color='brown')
#~ 
#~ plt.plot(goodk,[f_low]*len(goodk),zorder=3,color='brown',ls='solid')
#~ plt.plot(goodk,[-f_low]*len(goodk),zorder=3,color='brown',ls='solid')


vzcc_filt=np.real(np.fft.ifft2(vzkw))
vzcc_filt_norm=vzcc_filt/abs(vzcc_filt).max()

plotc.colorplot(vzcc_filt_norm,x=x,y=t,vmax=0.16,centerzero=True,subplot_index=next(g),ylim=[None,100],
xlabel="Distance (Mm)",ylabel="Time (min)",title="Unfiltered Time-Distance diagram for quiet Sun",usetex=True)

plotc.gca().text(-40,5.2,"Source",fontsize=15,horizontalalignment='center')

#~ plotc.drawline(x=2*15,color='brown')
#~ plotc.drawline(x=-2*15,color='brown')

plotc.gcf().set_size_inches(2*7.5,5)
plotc.savefig('/home/jishnu/project/magnetic_inversions/report/quiet_spectrum.eps')
#~ plt.show()

