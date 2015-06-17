from __future__ import division
import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.ndimage.filters import gaussian_filter1d
import os,fnmatch
import pyfits
import warnings

#######################################################################

def fitsread(f): 
    arr=pyfits.getdata(f)
    # If it is 2D, make it 3D. This adds an extra dimension at the end.
    # Bring it to the middle to match the general trend
    # Dimension will be (nz,ny,nx) after this
    if len(arr.shape)==2: arr=np.atleast_3d(arr).transpose(0,2,1)
    # Change dimension to (nx,ny,nz)
    return arr.transpose(2,1,0)

def fitswrite(f,arr): 
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        # Change dimensions from (nx,ny,nz) to (nz,ny,nx)
        arr=arr.transpose(2,1,0)
        pyfits.writeto(f,arr,clobber=True)

def get_iter_no():
    updatedir=os.path.join(datadir,"update")
    return len(fnmatch.filter(os.listdir(updatedir),'misfit_[0-9][0-9]'))
    
def check_BFGS_stencil(iterno,stencilBFGS):
    try: assert (iterno - 1 - stencilBFGS)>1
    except AssertionError:    
        print "BFGS stencil too wide, quitting."
        quit()

def get_number_of_sources():
    return len(np.loadtxt(os.path.join(datadir,'master.pixels')))

def filterx(kern):
    temp=np.zeros_like(kern)
    nx,ny,nz=kern.shape
    x=np.fft.fftfreq(nx)*nx
    smooth_x = 4
    x = np.exp(-x**2./(2.*smooth_x**2.))
    filtx=np.fft.rfft(x)
    filtx=np.atleast_3d(filtx).transpose(1,0,2)
    temp=np.fft.rfft(kern,axis=0)
    temp*=filtx
    temp=np.fft.irfft(temp,axis=0).real
    #temp=temp*nx*ny*nz # Remove if not comparing with fortran
    kern[:]=temp[:]

def filterz(arr,algo='spline',sp=1.0):
    nx,ny,nz=arr.shape
    z=np.arange(nz)
    b=arr.copy()
    
    if algo=='smooth':
        coeffs=np.zeros(6)
        coeffs[0]=0.75390625000000
        coeffs[1]=0.41015625000000
        coeffs[2]=-0.23437500000000
        coeffs[3]=0.08789062500000
        coeffs[4]=-0.01953125000000
        coeffs[5]=0.00195312500000
        sigmaz=4
        
        temp=np.zeros_like(arr)
        
        kst=1
        
        for k in xrange(5,nz-5):
            temp[:,:,k] =( coeffs[0]*arr[:,:,k] + 0.5*coeffs[1]*(arr[:,:,k-1] + arr[:,:,k+1]) 
                        +  0.5*coeffs[2]*(arr[:,:,k-2] + arr[:,:,k+2])  
                        +  0.5*coeffs[3]*(arr[:,:,k-3] + arr[:,:,k+3])  
                        +  0.5*coeffs[4]*(arr[:,:,k-4] + arr[:,:,k+4])  
                        +  0.5*coeffs[5]*(arr[:,:,k-5] + arr[:,:,k+5])  )
                        
        if (kst==1):

            temp[:,:,nz-5:nz-1] = 0.5 * (arr[:,:,nz-4:] + arr[:,:,nz-6:nz-2])
            temp[:,:,1:4] = 0.5 * (arr[:,:,:3] + arr[:,:,2:5])
            
        arr[:,:,kst:nz-kst+1]=temp[:,:,kst:nz-kst+1]
                
    elif algo=='spline':
        for x in xrange(nx):
            for y in xrange(ny):
                arrz=arr[x,y]
                arrzmax=arrz.max()
                arrz=arrz/arrzmax
                
                s=UnivariateSpline(z,arrz,s=sp)
                arr[x,y]=s(z)*arrzmax
                
    elif algo=='gaussian':
        arr[:]=gaussian_filter1d(arr,sigma=sp,axis=-1)
    
def antisymmetrize(arr):   arr[:]=0.5*(arr[:]-arr[::-1])

def twodigit(m): return str(m).zfill(2)
    
def iterminus(m): return twodigit(iterno-m)

def updatedir(filename): return os.path.join(datadir,'update',filename)

def rms(arr): return np.sqrt(np.sum(arr**2)/np.prod(arr.shape)) 

########################################################################



Rsun=695.9895 # Mm
    
datadir="/scratch/jishnu/magnetic/data"

steepest_descent = True
conjugate_gradient=True and not steepest_descent
LBFGS = True and not (steepest_descent or conjugate_gradient)

nx,ny,nz=256,1,300

itermax=25

iterno=get_iter_no()
back=np.loadtxt('polytrope')

num_src=get_number_of_sources()

array_shape=(nx,ny,nz)
totkern_c=np.zeros(array_shape)
totkern_psi=np.zeros(array_shape)
hess=np.zeros(array_shape)

for src in xrange(1,num_src+1):
    
    kern=fitsread(os.path.join(datadir,'kernel','kernel_c_'+twodigit(src)+'.fits'))
    totkern_c+=kern
    
    kern=fitsread(os.path.join(datadir,'kernel','kernel_vectorpsi_'+twodigit(src)+'.fits'))
    totkern_psi+=kern
    
    kern=fitsread(os.path.join(datadir,'kernel','hessian_'+twodigit(src)+'.fits'))
    hess+=abs(kern)

hess=hess*np.atleast_3d(back[:,2]).transpose(0,2,1)
hess = hess/abs(hess).max()
hess[hess<5e-3]=5e-3

kern = totkern_c/hess
filterx(kern)
filterz(kern,algo='smooth',sp=0.3)
totkern_c=kern


kern = totkern_psi/hess
filterx(kern)
antisymmetrize(kern)
filterz(kern,algo='gaussian',sp=2.0)
totkern_psi=kern

fitswrite(updatedir('gradient_c_'+iterminus(1)+'.fits'),totkern_c)
fitswrite(updatedir('gradient_vectorpsi_'+iterminus(1)+'.fits'),totkern_psi)

if (iterno==1) or steepest_descent:
    lastmodel_c=fitsread(updatedir('model_c_'+iterminus(1)+'.fits'))
    lastmodel_psi=fitsread(updatedir('model_vectorpsi_'+iterminus(1)+'.fits'))
    
    fitswrite(updatedir('update_c_'+iterminus(1)+'.fits'),totkern_c)
    fitswrite(updatedir('update_vectorpsi_'+iterminus(1)+'.fits'),totkern_psi)
    
    update_c = totkern_c 
    update_psi = totkern_psi
    
    if (iterno > 1): print 'Forcing steepest descent'

elif (iterno>1 and conjugate_gradient):
    lastmodel_c=fitsread(updatedir('model_c_'+iterminus(1)+'.fits'))
    lastmodel_psi=fitsread(updatedir('model_vectorpsi_'+iterminus(1)+'.fits'))
    
    grad_psi=fitsread(updatedir('gradient_vectorpsi_'+iterminus(1)+'.fits'))
    lastgrad_psi=fitsread(updatedir('gradient_vectorpsi_'+iterminus(2)+'.fits'))
    lastupdate_psi=fitsread(updatedir('update_vectorpsi_'+iterminus(2)+'.fits'))
    
    grad_c=fitsread(updatedir('gradient_c_'+iterminus(1)+'.fits'))
    lastgrad_c=fitsread(updatedir('gradient_c_'+iterminus(2)+'.fits'))
    lastupdate_c=fitsread(updatedir('update_c_'+iterminus(2)+'.fits'))
    
    con = np.sum(grad_psi*(grad_psi - lastgrad_psi))
    den = sum(lastgrad_psi**2.)
    con=con/den
    den=con
    if con<0: con=0
    print con,den
    
    update_c = totkern_c +  con* lastupdate_c
    update_psi = totkern_psi +  con* lastupdate_psi
    
    fitswrite(updatedir('update_c_'+iterminus(1)+'.fits'),update_c)
    fitswrite(updatedir('update_psi_'+iterminus(1)+'.fits'),update_psi)
    
elif (iterno>1 and LBFGS):
    
    stencilBFGS=4
    check_BFGS_stencil(iterno,stencilBFGS)
    minBFGS=iterno-1-stencilBFGS
    
    alph=np.zeros(iteration)
    
    model=fitsread(updatedir('model_c_'+iterminus(1)+'.fits'))
    grad_c=fitsread(updatedir('gradient_c_'+iterminus(1)+'.fits'))
    update_c=grad_c
    
    for i in xrange(iterno-2,minBFGS-1,-1):
        
        previterind=twodigit(i)
        lastmodel=fitsread(updatedir('model_c_'+previterind+'.fits'))
        lastgrad_c=fitsread(updatedir('gradient_c_'+previterind+'.fits'))
        
        alph[i] = np.sum((model-lastmodel)*update_c) /np.sum((model-lastmodel)*(grad_c-lastgrad_c))
        update_c = update_c - alph[i] * (grad_c - lastgrad_c) 
        
        model = lastmodel
        grad_c = lastgrad_c
        
        
    lastmodel=fitsread(updatedir('model_c_'+twodigit(minBFGS-1)+'.fits'))
    lastgrad_c=fitsread(updatedir('gradient_c_'+twodigit(minBFGS-1)+'.fits'))
    
    for i in xrange(minBFGS,iterno-1):
        model=fitsread(updatedir('model_c_'+twodigit(i)+'.fits'))
        grad_c=fitsread(updatedir('gradient_c_'+twodigit(i)+'.fits'))
        
        alph[i] = alph[i] - sum((grad_c-lastgrad_c)*update_c) /sum((model-lastmodel)*(grad_c-lastgrad_c))
        update_c = update_c + alph[i] * (model - lastmodel) 


        lastmodel = model
        lastgradc = gradc
        
psimax = abs(update_psi).max()

update_c = update_c/psimax 
update_psi = update_psi/psimax 

eps = 1e-4

psi_scale=rms(lastmodel_psi)

for i in xrange(1,6):
    update = lastmodel_c
    fitswrite(updatedir('test_c_'+str(i)+'.fits'), update)
    
    update = lastmodel_psi + eps * i * psi_scale * update_psi
    fitswrite(updatedir('test_vectorpsi_'+str(i)+'.fits'), update)
