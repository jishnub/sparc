from __future__ import division
import numpy as np
import plotc
import pyfits as pf



for source in xrange(1,9):
	master_pixels=np.loadtxt('/home/jishnu/project/magnetic_inversions/working/data/master.pixels').astype(int)
	src_loc=master_pixels[source-1]
	source=str(source)
	label=source+" at "+str(src_loc)
	source=source.zfill(2)	

	vzcc=np.squeeze(pf.getdata('/home/jishnu/project/magnetic_inversions/working/data/forward'+source+'/vz_cc.fits'))
	quiet=np.squeeze(pf.getdata('/home/jishnu/project/magnetic_inversions/working/quiet.fits'))

	nxt=vzcc.shape
	nx=nxt[1]
	nt=nxt[0]

	fmode=np.zeros((nx,1,nt))
	xlength=400*10**8
	L=xlength*10.**(-8.)*nx/(nx-1.)
	dx=L/nx
	dt = 30
	pi=np.pi
	cos=np.cos

	active_source_index=abs(vzcc[17]).argmax()
	quiet_source_index=abs(quiet[17]).argmax()
	#~ print "Active source at",active_source_index
	#~ print "Quiet source at",quiet_source_index

	kf = np.fft.fftfreq(nx,dx)*2*pi
	ff = np.fft.fftfreq(nt,dt)

	k = abs(kf)
	w = abs(ff) * 2.*pi
	f = w/(2.*pi)*1e3
	ff=ff*1e3

	kax,fax=np.meshgrid(k,f)

	x=np.linspace(-xlength/2,xlength/2,nx)/1e8
	t=np.linspace(0,nt*dt/60,nt)
	Xax,Tax=np.meshgrid(x,t)

	#~ Define f mode filter
	Poly=np.zeros(3)
	Poly[0]=0.0010
	Poly[1]=0.0030
	Poly[2]=-0.0006
	f_mode_const=0.00293795
	f0=0.7*f_mode_const*k**0.5*1e3
	f1=(Poly[0] + Poly[1]*k +Poly[2]*k**2.)*1e3



	for i in xrange(nx):
	  delta = (f1[i] - f0[i])
	  for j in xrange(nt):
		d = f[j] - f0[i]  
		if d>0 and d < delta:
		  fmode[i,0,j] = 0.5*(1.+cos(pi*(d-delta*0.5)/(delta*0.5)))
	df = 0.5
	f_low = 1.1
	for j in xrange(nt):
	  if (f[j] < f_low+df):
		fmode[:,0,j] = fmode[:,0,j]*0.5*(1.+cos(pi*(f[j]-(f_low+df))/df) )
	  if (f[j] < f_low): fmode[:,0,j] = 0.

	fmode=np.squeeze(fmode).T

	#~ Define p mode filter
	Poly[0]=0.0011
	Poly[1]=0.0052
	Poly[2]=-0.0013
	f_mode_const=0.0033
	pmode=np.zeros((nx,1,nt))
	f0=1.2*f_mode_const*k**0.5*1e3
	f1=(Poly[0] + Poly[1]*k +Poly[2]*k**2.)*1e3
	f_low = 1.6

	for i in xrange(nx):
	  delta = (f1[i] - f0[i])
	  for j in xrange(nt):
		d = f[j] - f0[i]
		if d< delta and d>0:
		  pmode[i,0,j] = 0.5*(1.+cos(pi*(abs(d)-abs(delta)*0.5)/(abs(delta)*0.5)))

	for j in xrange(nt):
	  if (f[j] <f_low): pmode[:,0,j] = 0.
	  if (f[j] < f_low+df):
		pmode[:,0,j] = pmode[:,0,j]*0.5*(1.+cos(pi*(f[j]-(f_low+df))/df) )

	pmode=np.squeeze(pmode).T

	filtcode="f"
	if filtcode=="f":
		filt=fmode
	elif filtcode=="p":
		filt=pmode

	vzcc_filt=np.real(np.fft.ifft2(np.fft.fft2(vzcc)*filt))
	vzcc_filt_norm=vzcc_filt/abs(vzcc_filt).max()

	quiet_filt=np.real(np.fft.ifft2(np.fft.fft2(quiet)*filt))
	quiet_filt_norm=quiet_filt/abs(quiet_filt).max()

	#~ Active
	pmode_vel_a=1.465
	fmode_vel_a=1.3

	if filtcode=="f":
		mode_vel_a=fmode_vel_a
	elif filtcode=="p":
		mode_vel_a=pmode_vel_a

	#~ Quiet
	pmode_vel_q=1.465
	fmode_vel_q=0.8

	if filtcode=="f":
		mode_vel_q=fmode_vel_q
	elif filtcode=="p":
		mode_vel_q=pmode_vel_q

	dist_low=-50
	dist_high=50
	dist_range=np.arange(dist_low,dist_high)
	rms_ratio=[]
	max_ratio=[]
	
	for dist_ind in dist_range:
		
		#~ plotc.figure()
		#~ g=plotc.gridlist(1,2)
		
		#~ plotc.colorplot(quiet_filt_norm,y=t,subplot_index=next(g))
		
		low_t_q=abs(x[quiet_source_index+dist_ind]-x[quiet_source_index])/mode_vel_q
		low_t_q_ind=np.where(t>low_t_q)[0][0]
		window_q=40
		high_t_q=abs(x[quiet_source_index+dist_ind]-x[quiet_source_index])/mode_vel_q+window_q
		high_t_q_ind=low_t_q_ind+window_q*60//dt
		
		#~ print low_t_q,low_t_q_ind,high_t_q,high_t_q_ind
		
		#~ plotc.drawline(y=low_t_q,color='g')
		#~ plotc.drawline(y=high_t_q,color='g')
		#~ plotc.drawline(x=dist_ind+quiet_source_index,color='g')
		
		#~ plotc.colorplot(vzcc_filt_norm,y=t,subplot_index=next(g),vmin=-0.1,centerzero=True)
		
		low_t_a=abs(x[active_source_index+dist_ind]-x[active_source_index])/mode_vel_a
		low_t_a_ind=np.where(t>low_t_a)[0][0]
		window_a=40
		high_t_a=abs(x[active_source_index+dist_ind]-x[active_source_index])/mode_vel_a+window_a
		high_t_a_ind=low_t_a_ind+window_a*60//dt
		
		#~ plotc.drawline(y=low_t_a,color='b')
		#~ plotc.drawline(y=high_t_a,color='b')
		#~ plotc.drawline(x=dist_ind+active_source_index,color='b')

		#~ plotc.figure()
		#~ plotc.plot1D(vzcc_filt[:,active_source_index+dist_ind],x=t,label="active")
		#~ plotc.plot1D(quiet_filt[:,quiet_source_index+dist_ind],x=t,label="quiet")
	#~ 
		#~ plotc.drawline(x=low_t_a,color='b')
		#~ plotc.drawline(x=high_t_a,color='b')
		#~ plotc.drawline(x=low_t_q,color='g')
		#~ plotc.drawline(x=high_t_q,color='g')
		#~ plotc.legend()
		#~ plotc.show()

		act_data=vzcc_filt[low_t_a_ind:high_t_a_ind,dist_ind+active_source_index]
		quiet_data=quiet_filt[low_t_q_ind:high_t_q_ind,dist_ind+quiet_source_index]

		act_max=abs(act_data).max()
		quiet_max=abs(quiet_data).max()
		
		act_rms=np.sqrt(sum(act_data**2)/len(act_data))
		quiet_rms=np.sqrt(sum(quiet_data**2)/len(quiet_data))

		#~ print quiet_max,act_max,quiet_max/act_max
		max_ratio.append(quiet_max/act_max)
		rms_ratio.append(quiet_rms/act_rms)

		#~ plotc.figure()
		#~ plotc.plot1D(act_data)
		#~ plotc.plot1D(quiet_data)
		#~ plotc.show()

	plotc.plot1D(rms_ratio,x=dist_range,label=label,subplot_index="121",title="RMS ratio")
	plotc.plot1D(max_ratio,x=dist_range,label=label,subplot_index="122",title="Max ratio")
	
#~ plotc.legend()
plotc.show()
