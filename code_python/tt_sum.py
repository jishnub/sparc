import numpy as np
import matplotlib.pyplot as plt

tt_file='/home/jishnu/project/magnetic_inversions/working/data/update/linesearch_00'

tdata=np.loadtxt(tt_file)

ones=np.where(tdata[:,0]==1)[0]
print ones

tt=np.zeros(len(ones))
for i in xrange(len(ones)-1):
	tt[i]=sum(tdata[ones[i]:ones[i+1],2])
	#~ print "Ignoring",int(tdata[ones[i]+7,0])
	#~ tt[i]=tt[i]-tdata[ones[i]+7,2]
	#~ if i==0: tt[i]=tt[i]-tdata[ones[i],2]-tdata[ones[i]+3,2]
	#~ if i==2: tt[i]=tt[i]-tdata[ones[i]+3,2]
	print int(tdata[ones[i],0]),':',int(tdata[ones[i+1]-1,0]),"\t",tt[i]

#~ print "Ignoring",int(tdata[ones[-1]+7,0])
#~ tt[-1]=sum(tdata[ones[-1]:,2])-tdata[ones[-1]+7,2]
tt[-1]=sum(tdata[ones[-1]:,2])
print int(tdata[ones[-1],0]),':',int(tdata[-1,0]),"\t",tt[-1]
