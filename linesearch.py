import os,sys,shutil,glob,time,re,subprocess

env=dict(os.environ, MPI_TYPE_MAX="1280280")

codedir="/home/jishnu/sparc/"
data="/scratch/jishnu/magnetic/data"

procno=int(os.environ["PBS_VNODENUM"])
with open(os.path.join(data,'master.pixels'),'r') as mp:
    nmasterpixels=sum(1 for _ in mp)

total_no_of_linesearches=5
total_no_of_jobs=nmasterpixels*total_no_of_linesearches

if procno>=total_no_of_jobs: quit()

linesearch_no=procno/nmasterpixels+1
src_no=procno%nmasterpixels+1

linesearch_no=str(linesearch_no).zfill(2)
src=str(src_no).zfill(2)

def compute_forward(linesearch_no,src):
    
    Spectral=os.path.join(codedir,"Spectral")
    Instruction=os.path.join(codedir,"Instruction_"+src+"_"+linesearch_no)
    forward = os.path.join(data,"forward_src"+src+"_ls"+linesearch_no)
    
    shutil.copyfile(Spectral,Instruction)
    
    sparccmd="/home/apps/openmpi-1.6.5/bin/mpiexec -np 1 ./sparc "+src+" "+linesearch_no
    with open(os.path.join(data,forward,"out_linesearch_"+linesearch_no),'w') as outfile:
        fwd=subprocess.call(sparccmd.split(),stdout=outfile,env=env,cwd=codedir)
    
    shutil.copyfile(os.path.join(forward,"vz_cc.fits"),os.path.join(forward,"vz_cc_"+str(linesearch_no)+".fits"))
    
    partialfiles=glob.glob(os.path.join(forward,"*partial*"))
    for f in partialfiles: os.remove(f)
    
    fullfiles=glob.glob(os.path.join(forward,"*full*"))
    for f in fullfiles: os.remove(f)
    
    return 0
    
compute_forward(linesearch_no,src)
