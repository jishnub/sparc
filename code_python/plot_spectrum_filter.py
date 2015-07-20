from __future__ import division
import numpy as np
import plotc
import pyfits
import os

Rsun=695.8 # Mm

username=os.environ['USER']
homedir=os.environ['HOME']
datafile=os.path.join("/scratch",username,"magnetic","data","data","01.fits")
data=np.squeeze(pyfits.getdata(datafile))

codedir=os.path.join(homedir,"sparc")
fmode_filter_file=os.path.join(codedir,"fmode_filter.fits")
pmode_filter_file=os.path.join(codedir,"pmode_filter.fits")

fmode_filter=np.squeeze(pyfits.getdata(fmode_filter_file))
pmode_filter=np.squeeze(pyfits.getdata(pmode_filter_file))

nt,nx=fmode_filter.shape

data_fft=np.fft.fft2(data)

data_ampspec=abs(data_fft)

data_f_filtered=data_fft*fmode_filter
data_p_filtered=data_fft*pmode_filter

data_f_filtered_ampspec=abs(data_f_filtered)
data_p_filtered_ampspec=abs(data_p_filtered)


x=np.linspace(-200,200,nx,endpoint=False) # Mm
dx=x[1]-x[0]
k=np.fft.fftfreq(nx,dx)*2*np.pi # Mm^-1

dt=30.0 # Seconds
nu=np.fft.fftfreq(nt,dt)*1e3 # MHz

# f-mode

f_mode_const=2.93795*0.7
f0=f_mode_const*np.sqrt(abs(k))

Poly=np.zeros(3)
Poly[0]=1.1
Poly[1]=5.2
Poly[2]=-1.3
f1=Poly[0] + Poly[1]*k +Poly[2]*k**2.

f_low = 1.1 # MHz

plotc.plt.figure()
plotc.spectrumplot(data_ampspec,x=k*Rsun,y=nu,sp=121,xr=[0,None],yr=[0,6],
                    axes_properties=dict(title="Full spectrum",xlabel=r"$k R_{\odot}$",ylabel=r"$ \nu $ (mHz)"))

plotc.plt.plot(k[:len(k)//2]*Rsun,f0[:len(k)//2],color="k")
plotc.plt.plot(k[:len(k)//2]*Rsun,f1[:len(k)//2],color='k')
plotc.plt.plot(k[:len(k)//2]*Rsun,[f_low]*len(k[:len(k)//2]))
plotc.spectrumplot(data_f_filtered_ampspec,x=k*Rsun,y=nu,sp=122,xr=[0,None],yr=[0,6],
                    axes_properties=dict(title="f-mode",xlabel=r"$k R_{\odot}$",ylabel=r"$ \nu $ (mHz)"))
                    
# p-mode

f_mode_const=1.2*3.3

Poly[0]=1.1
Poly[1]=5.2
Poly[2]=-1.3

f0=f_mode_const*abs(k)**0.5
f1=Poly[0] + Poly[1]*k +Poly[2]*k**2.
f_low = 1.6 # mHz

plotc.plt.figure()
plotc.spectrumplot(data_ampspec,x=k*Rsun,y=nu,sp=121,xr=[0,None],yr=[0,6],
                    axes_properties=dict(title="Full spectrum",xlabel=r"$k R_{\odot}$",ylabel=r"$ \nu $ (mHz)"))


plotc.plt.plot(k[:len(k)//2]*Rsun,f0[:len(k)//2],color="k")
plotc.plt.plot(k[:len(k)//2]*Rsun,f1[:len(k)//2],color='k')
plotc.plt.plot(k[:len(k)//2]*Rsun,[f_low]*len(k[:len(k)//2]))
plotc.spectrumplot(data_p_filtered_ampspec,x=k*Rsun,y=nu,sp=122,xr=[0,None],yr=[0,6],
                    axes_properties=dict(title="p-mode",xlabel=r"$k R_{\odot}$",ylabel=r"$ \nu $ (mHz)"))
                    
plotc.plt.show()

