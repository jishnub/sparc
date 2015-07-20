import os,shutil,glob,re,subprocess

env=dict(os.environ, MPI_TYPE_MAX="1280280")

codedir=os.path.dirname(os.path.abspath(__file__))
HOME=os.environ["HOME"]

configvars={}
with open(os.path.join(codedir,"varlist.sh")) as myfile:
    for line in myfile:
        name,var=line.partition("=")[::2]
        configvars[name.strip()]=var.strip().strip('"')
        
data=configvars['directory']

iterno=len([f for f in os.listdir(os.path.join(data,'update')) if re.match(r'misfit_[0-9]{2}$',f)])
itername=str(iterno).zfill(2)

with open(os.path.join(data,'master.pixels'),'r') as mp:
    nmasterpixels=sum(1 for _ in mp)

procno=int(os.environ["PBS_VNODENUM"])

if procno>=nmasterpixels: quit()

src=str(procno+1).zfill(2)

def copyfile(pathA,pathB):
    try:
        shutil.copyfile(pathA,pathB)
    except IOError:
        print "Error copying",pathA        

def compute_forward_adjoint_kernel(src):

    forward="forward_src"+src+"_ls00"
    adjoint="adjoint_src"+src
    kernel="kernel"

    ttname="vz_cc_src"+src+".fits"
    tdiff0name="ttdiff_src"+src+".fmode"
    tdiff1name="ttdiff_src"+src+".p1mode"
    Spectral=os.path.join(codedir,"Spectral")
    Adjoint=os.path.join(codedir,"Adjoint")
    Instruction=os.path.join(codedir,"Instruction_src"+src+"_ls00")
    
    
    mpipath=os.path.join(HOME,"anaconda/bin/mpiexec")
    sparccmd=mpipath+" -np 1 ./sparc "+src+" 00"
    
    def inforward(filename): return os.path.join(data,forward,filename)
    def inadjoint(filename): return os.path.join(data,adjoint,filename)
    def inkernel(filename): return os.path.join(data,kernel,filename)
    def indata(filename): return os.path.join(data,filename)
    def instatus(filename): return os.path.join(data,"status",filename)
    def inttiter(filename): return os.path.join(data,"tt","iter"+itername,filename)
    
    ####################################################################
    #~ Forward
    ####################################################################
    
    if not os.path.exists(instatus(forward)):

        copyfile(Spectral,Instruction)

        with open(inforward("out_"+forward),'w') as outfile:
            fwd=subprocess.call(sparccmd.split(),stdout=outfile,env=env,cwd=codedir)
        
        if not os.path.exists(os.path.join(data,"tt","iter"+itername)):
            os.makedirs(os.path.join(data,"tt","iter"+itername))
                        
        copyfile(inforward("vz_cc.fits"),inforward("vz_cc_00.fits"))
        
        copyfile(inforward("vz_cc.fits"),inttiter(ttname))
                        
        copyfile(inforward("ttdiff.0"),inttiter(tdiff0name))
                        
        copyfile(inforward("ttdiff.1"),inttiter(tdiff1name))
    
    ####################################################################
    #~ Adjoint
    ####################################################################

    if not os.path.exists(instatus(adjoint)):
        copyfile(Adjoint,Instruction)

        with open(inadjoint("out_"+adjoint),'w') as outfile:
            adj=subprocess.call(sparccmd.split(),stdout=outfile,env=env,cwd=codedir)
            
    ####################################################################
    #~ Kernel
    ####################################################################

    if not os.path.exists(instatus(kernel+src)):

        with open(inkernel("out_kernel"+src),'w') as outfile:
            kern=subprocess.call(sparccmd.split(),stdout=outfile,env=env,cwd=codedir)
    
    return 0
    

compute_forward_adjoint_kernel(src)
