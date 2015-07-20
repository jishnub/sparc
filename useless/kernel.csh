#PBS -S /bin/csh

# Job Name:
#PBS -N COMPUTE_KERNELS
#PBS -q debug

#PBS -l select=1:ncpus=1
#PBS -l walltime=00:15:00

module add mpi-sgi/mpt.1.26
setenv MPI_TYPE_MAX 1280280

mpiexec -np 1 ./sparc > outputkern
rm core.*

