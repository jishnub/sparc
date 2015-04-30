#! /bin/csh

touch running_full
setenv MPI_TYPE_MAX 1280280
set directory = /home/jishnu/project/magnetic_inversions/working3/data
set nmasterpixels = `wc -l < $directory/master.pixels`

find $directory -name "linesearch" -exec rm -f {} \; 
find $directory -name "compute_data" -exec rm -f {} \; 
find $directory -name "compute_synth" -exec rm -f {} \; 

set iter=`find $directory/update -name 'misfit_[0-9][0-9]'|wc -l`
set itername=`printf "%02d" $iter`

foreach i (`seq -f "%02g" 1 $nmasterpixels`)
  
  set forwardname = forward$i
  set adjointname = adjoint$i
  set kernelname = kernel$i
  set ttname = vz_cc_src$i.fits
  set tdiff0name = ttdiff_src$i.0
  set tdiff1name = ttdiff_src$i.1

  if (! -e $directory/status/$forwardname) then
   echo "Running "$forwardname
   cp Spectral Instruction
   mpiexec -np 1 ./sparc > $directory/$forwardname/out$forwardname
   mkdir -p $directory/tt/iter$itername
    cp $directory/$forwardname/vz_cc.fits $directory/tt/iter$itername/$ttname
    cp $directory/$forwardname/ttdiff.0 $directory/tt/iter$itername/$tdiff0name
    cp $directory/$forwardname/ttdiff.1 $directory/tt/iter$itername/$tdiff1name
  endif

  if (! -e $directory/status/$adjointname) then
   echo "Running "$adjointname
   cp Adjoint Instruction
   mpiexec -np 1 ./sparc > $directory/$adjointname/out$adjointname
  endif

  if (! -e $directory/status/$kernelname) then
   echo "Running "$kernelname
   mpiexec -np 1 ./sparc > $directory/$adjointname/out$kernelname
  endif


  find . -name "core.*" -exec rm -f {} \; 
  find . -name "fort.*" -exec rm -f {} \; 

end
  
find $directory/status -name "forward*" -exec rm -f {} \;
find $directory/status -name "adjoint*" -exec rm -f {} \;
find $directory/status -name "kernel*" -exec rm -f {} \;

cp $directory/kernel/misfit $directory/update/misfit_$itername
cp $directory/kernel/misfit_all $directory/update/misfit_all_$itername
#cp vx.fits $directory/update/vx_$itername.fits
#cp vz.fits $directory/update/vz_$itername.fits
find $directory/kernel -name "misfit" -exec rm -f {} \;
find $directory/kernel -name "misfit_all" -exec rm -f {} \;
find . -name "fort.*" -exec rm -f {} \; 

rm running_full
