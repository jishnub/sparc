#PBS -S /bin/csh

# Job Name:
#PBS -N COMPUTE_GRAD
#PBS -q devel

#PBS -l select=1:ncpus=1:model=san
#PBS -l walltime=02:00:00

module add mpi-sgi/mpt.1.26
setenv MPI_TYPE_MAX 1280280
set directory = /nobackupp2/shanasog/classic
set nmasterpixels = `wc -l < $directory/master.pixels`
@ i = 1

find $directory -name "linesearch" -exec rm -f {} \; 

@ iter = -1
@ log = 1

while ($log == 1)
 
 @ iter++
 if ($iter < 10) then
  set itername = 0$iter
 endif

 if ($iter >= 10) then
  set itername = $iter
 endif

 @ log = 0
 if (-e $directory/update/misfit_$itername) then
  @ log = 1 
 endif

end

 if (! -d $directory/tt/$itername) then
  mkdir $directory/tt/$itername
 endif

while ($i <= $nmasterpixels) 
 
 if ($i < 10) then
  set forwardname = forward0$i
  set adjointname = adjoint0$i
  set kernelname = kernel0$i
  set ttname = 0$i.fits
 else
  set forwardname = forward$i
  set adjointname = adjoint$i
  set kernelname = kernel$i
  set ttname = $i.fits
 endif

  if (! -e $directory/status/$forwardname) then
   cp Spectral Instruction
   mpiexec -np 1 ./sparc > $directory/$forwardname/out$forwardname
   cp $directory/$forwardname/vz_cc.fits $directory/tt/$itername/$ttname
  endif

  if (! -e $directory/status/$adjointname) then
   cp Adjoint Instruction
   mpiexec -np 1 ./sparc > $directory/$adjointname/out$adjointname
  endif

  if (! -e $directory/status/$kernelname) then
   mpiexec -np 1 ./sparc > $directory/$adjointname/out$kernelname
  endif


  find . -name "core.*" -exec rm -f {} \; 
  find . -name "fort.*" -exec rm -f {} \; 

  if ($i == $nmasterpixels) then
   find $directory/status -name "forward*" -exec rm -f {} \;
   find $directory/status -name "adjoint*" -exec rm -f {} \;
   find $directory/status -name "kernel*" -exec rm -f {} \;
  endif

  @ i++
end
  cp $directory/kernel/misfit $directory/update/misfit_$itername
  find $directory/kernel -name "misfit" -exec rm -f {} \;
