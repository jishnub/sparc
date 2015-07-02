import numpy as np
import plotc

l1=np.loadtxt('/home/jishnu/project/magnetic_inversions/working_2/data/update/linesearch_00')
l2=np.loadtxt('/home/jishnu/project/magnetic_inversions/working_4/data/update/linesearch_00')

line1=[sum(l1[8*i:8*(i+1),2]) for i in xrange(0,5)]
line2=[sum(l2[8*i:8*(i+1),2]) for i in xrange(0,5)]

line1.extend(line2)
line_full=np.array(line1)

epsilon=np.linspace(1e-4,10e-4,10)

plotc.plot1D(line_full,x=epsilon,yr=[6.1,7],xr=[0,1.1e-3],marker='o',color='g',
xlabel="Step size $\epsilon$",ylabel="Misfit",title="First iteration: line search",usetex=True)
plotc.figuresize(6.5,5)
plotc.gca().annotate("Optimal\nStep Size", xy=(0.7e-3, 6.22), xytext=(0.7e-3, 6.6),
    arrowprops=dict(arrowstyle="->"),fontsize=15,horizontalalignment='center')
#~ plotc.show()
plotc.savefig('/home/jishnu/project/magnetic_inversions/report/misfit_iter1.eps')
