import os,sys,shutil,glob,subprocess


env=dict(os.environ, MPI_TYPE_MAX="1280280")
codedir=os.path.dirname(os.path.abspath(__file__))

configvars={}
with open(os.path.join(codedir,"varlist.sh")) as myfile:
    for line in myfile:
        name,var=line.partition("=")[::2]
        configvars[name.strip()]=var.strip().strip('"')

data=configvars['directory']

procno=int(os.environ["PBS_VNODENUM"])
HOME=os.environ["HOME"]

with open(os.path.join(data,'master.pixels'),'r') as mp:
    nmasterpixels=sum(1 for _ in mp)

if procno>=nmasterpixels: quit()

src=str(procno+1).zfill(2)

def copyfile(pathA,pathB):
    try:
        shutil.copyfile(pathA,pathB)
    except IOError:
        print "Error copying",pathA,"to",pathB
        
        
def movefile(pathA,pathB):
    try:
        shutil.move(pathA,pathB)
    except IOError:
        print "Error moving",pathA,"to",pathB

def compute_data(src):

    forward="forward_src"+src+"_ls00"
    Spectral=os.path.join(codedir,"Spectral")
    Instruction=os.path.join(codedir,"Instruction_src"+src+"_ls00")
    
    def inforward(filename): return os.path.join(data,forward,filename)
    
    copyfile(Spectral,Instruction)
    
    mpipath=os.path.join(HOME,"anaconda/bin/mpiexec")
    sparccmd=mpipath+" -np 1 ./sparc "+src+" 00"
    
    with open(inforward("out_data_forward"),'w') as outfile:
        fwd=subprocess.call(sparccmd.split(),stdout=outfile,env=env,cwd=codedir)

    partialfiles=glob.glob(inforward("*partial*"))
    for f in partialfiles: os.remove(f)
    
    fullfiles=glob.glob(inforward("*full*"))
    for f in fullfiles: os.remove(f)
    
    if not os.path.exists(os.path.join(data,"tt","data")):
        os.makedirs(os.path.join(data,"tt","data"))
    
    if os.path.exists(inforward("vz_cc.fits")):
        movefile(inforward("vz_cc.fits"),inforward("data.fits"))
        copyfile(inforward("data.fits"),os.path.join(data,"tt","data","data"+src+".fits"))
        copyfile(inforward("data.fits"),os.path.join(data,"data",src+".fits"))
    else:
        print "vz_cc.fits not found. Check the forward calculation"

compute_data(src)

file_to_remove=os.path.join(data,"status","forward_src"+src+"_ls00")
if os.path.exists(file_to_remove): os.remove(file_to_remove)

