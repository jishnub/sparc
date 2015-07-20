from __future__ import division
import numpy as np
import pyfits as pf
import plotc
from mpldatacursor import datacursor
from scipy.optimize import curve_fit
import os

Lx=400 #Mm
dt=30/60 #minutes
username=os.environ['USER']
codedir=os.path.join('/home',username,'sparc')

configvars={}
with open(os.path.join(codedir,"varlist.sh")) as myfile:
    for line in myfile:
        name,var=line.partition("=")[::2]
        configvars[name.strip()]=var.strip().strip('"')

datadir=configvars['directory']

master_pixels=np.loadtxt(os.path.join(datadir,'master.pixels'))
Nsources=len(master_pixels)

wavespeed=np.zeros((Nsources,2))

p1mode=True
fmode=False

for sourceno in xrange(Nsources):
    sourceloc=master_pixels[sourceno]
    sourcedir='forward_src'+str(sourceno+1).zfill(2)+'_ls00'

    vzcc=np.squeeze(pf.getdata(os.path.join(datadir,sourcedir,'data.fits')))
    Nx=vzcc.shape[1]
    sourcepix=int(np.floor(sourceloc/Lx*Nx)+Nx/2-1)
    if fmode:
        
        
        plotc.plt.figure()

        fmode_filter=np.squeeze(pf.getdata(os.path.join(codedir,'fmode_filter.fits')))

        vzcc_f=np.real(np.fft.ifft2(np.fft.fft2(vzcc)*fmode_filter))
        vzcc_f/=abs(vzcc_f).max()
        plotc.colorplot(vzcc_f,vmax=0.01,centerzero=True)
        datacursor(display='multiple',draggable='True')



        Npoints=5
        xpix=np.array([-30-ind*6 for ind in xrange(Npoints)])+sourcepix
        xpoints=(xpix-127)/256*400

        plotc.draw_vlines(xpix)

        t_cutoff=np.zeros((Npoints,2))

        #~ Line fit to bottom of wavepacket
        def bot_time(x):
            return (89.2-161)/(89.2-60.1)*(x-sourcepix)+20

        #~ Line fit to top of wavepacket
        def top_time(x):
            return (258-149)/(68.5-94.3)*(x-sourcepix)+60

        t_cutoff[:,0]=bot_time(xpix)
        t_cutoff[:,1]=top_time(xpix)
        
        plotc.plt.plot(xpix,t_cutoff[:,0],'b-')
        plotc.plt.plot(xpix,t_cutoff[:,1],'g-')

        xdat=vzcc_f[:,xpix].T

        def local_max_index(arr):
            return (np.diff(np.sign(np.diff(arr))) < 0).nonzero()[0] + 1

        peakcenter=[]
        peakamplitude=[]
        peakwidth=[]

        def gaus(x,a,x0,sigma):
            return a*np.exp(-(x-x0)**2/(2*sigma**2))


        for xind,tseries in enumerate(xdat):
            tp=local_max_index(tseries)
            tp=tp[np.where(np.logical_and(tp<t_cutoff[xind,1],tp>t_cutoff[xind,0]))[0]]
            yp=tseries[tp]    
            
            plotc.plt.plot([xpix[xind]]*len(tp),tp,'o',color='orange')

            popt,_ = curve_fit(gaus,tp,yp,p0=[yp[len(yp)//2],tp[len(tp)//2],(tp[-1]-tp[0])*0.4])
            
            peakamplitude.append(popt[0])
            peakcenter.append(popt[1])
            peakwidth.append(popt[2])
            
        peakamplitude=np.array(peakamplitude)
        peakcenter=np.array(peakcenter)*dt
        peakwidth=np.array(peakwidth)*dt
            
        v,_=np.polyfit(peakcenter,xpoints,1)

        print "Source",sourceno+1,"velocity of f mode",abs(v),"Mm/min"
        
        wavespeed[sourceno,0]=abs(v)
        
        plotc.plt.show()

    ####################################################################
    #~ p-mode
    
    if p1mode:
        plotc.plt.figure()
        
        pmode_filter=np.squeeze(pf.getdata(os.path.join(codedir,'pmode_filter.fits')))

        vzcc_p=np.real(np.fft.ifft2(np.fft.fft2(vzcc)*pmode_filter))
        vzcc_p/=abs(vzcc_p).max()
        plotc.colorplot(vzcc_p,vmax=0.1,centerzero=True)
        datacursor(display='multiple',draggable='True')



        Npoints=5
        xpix=np.array([-30-ind*6 for ind in xrange(Npoints)])+sourcepix
        xpoints=(xpix-127)/256*400

        plotc.draw_vlines(xpix)

        t_cutoff=np.zeros((Npoints,2))

        #~ Line fits are in pixels
        #~ Line fit to bottom of wavepacket
        #~ def bot_time(x):
            #~ return (19.7-190)/(98.1-16.1)*(x-sourcepix)
        def bot_time(x):
            return (101-61)/(73.7-92.4)*(x-92.4)+61

        #~ Line fit to top of wavepacket
        def top_time(x):
            return (153-258)/(70.3-28)*(x-sourcepix)+58

        t_cutoff[:,0]=bot_time(xpix)
        t_cutoff[:,1]=top_time(xpix)
        
        plotc.plot(xpix,t_cutoff[:,0],'b-')
        plotc.plot(xpix,t_cutoff[:,1],'g-')


        xdat=vzcc_p[:,xpix].T

        def local_max_index(arr):
            return (np.diff(np.sign(np.diff(arr))) < 0).nonzero()[0] + 1

        peakcenter=[]
        peakamplitude=[]
        peakwidth=[]
        failind=[]

        def gaus(x,a,x0,sigma):
            return a*np.exp(-(x-x0)**2/(2*sigma**2))


        for xind,tseries in enumerate(xdat):
            tp=local_max_index(tseries)
            tp=tp[np.where(np.logical_and(tp<t_cutoff[xind,1],tp>t_cutoff[xind,0]))[0]]
            yp=tseries[tp]    

            plotc.plot([xpix[xind]]*len(tp),tp,'o',color='orange')
            
            try:
                popt,_ = curve_fit(gaus,tp,yp,p0=[yp[len(yp)//2],tp[len(tp)//2],(tp[-1]-tp[0])*0.4])
            except RuntimeError:
                print "Gaussian fit failed for x",xpix[xind]
                failind.append(xind)
                continue
            
            peakamplitude.append(popt[0])
            peakcenter.append(popt[1])
            peakwidth.append(popt[2])
            
        peakamplitude=np.array(peakamplitude)
        peakcenter=np.array(peakcenter)*dt
        peakwidth=np.array(peakwidth)*dt
        
        try:
            v,_=np.polyfit(peakcenter,np.delete(xpoints,failind),1)
            print "Source",sourceno+1,"velocity of p mode",abs(v),"Mm/min"
            wavespeed[sourceno,1]=abs(v)
        except RuntimeError:
            print "Fit did not converge for source",sourceno
            continue
        
        plotc.show()

#~ np.savetxt(os.path.join(datadir,'wavespeed'),wavespeed,fmt="%10.8f")



