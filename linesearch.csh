#! /bin/csh

#~ module add mpi-sgi/mpt.1.26
#~ setenv MPI_TYPE_MAX 1280280

touch running_linesearch
set directory = /home/jishnu/project/magnetic_inversions/working3/data
set nmasterpixels = `wc -l < "$directory/master.pixels"`

set nmas = `printf "%02d" $nmasterpixels`

@ j = 1

set iter=`find $directory/update -name 'linesearch_[0-9][0-9]'|wc -l`
set itername=`printf "%02d" $iter`

##echo "Make sure misfit has been moved"

touch $directory/linesearch
echo "Starting iterations"
while ($j <= 5)
 echo $j
 if (! -e $directory/update/tested$j) then

  if (-e $directory/status/forward$nmas) then
   find $directory/status -name "forward*" -exec rm -f {} \; 
  endif

  
  foreach src (`seq -f "%02g" 1 $nmasterpixels`) 
    echo $src
    set forwardname = $directory/forward$src

    cp Spectral Instruction
    cp $directory/update/test_c_$j.fits $directory/model_c.fits
    cp $directory/update/test_vectorpsi_$j.fits $directory/model_vectorpsi.fits
    mpiexec -np 1 ./sparc > $forwardname/output
    cp $forwardname/vz_cc.fits $forwardname/vz_cc_$j.fits
    find $forwardname -name "*partial*" -exec rm -f {} \; 
    find $forwardname -name "*full*" -exec rm -f {} \; 
  end
  
  touch $directory/update/tested$j
 endif

 @ j++

end

if (-e $directory/update/tested5) then
 find $directory/update -name "tested*" -exec rm -f {} \; 
 #find $directory -name "model_c.fits" -exec rm -f {} \; 
 #find $directory -name "model_psi.fits" -exec rm -f {} \; 
 find $directory -name "update.fits" -exec rm -f {} \; 
 find $directory -name "linesearch" -exec rm -f {} \; 
 find $directory/status -name "forward*" -exec rm -f {} \;
 cp $directory/kernel/misfit $directory/update/linesearch_$itername
 cp $directory/kernel/misfit_all $directory/update/linesearch_all_$itername
 find $directory/kernel -name "misfit" -exec rm -f {} \;
 find $directory/kernel -name "misfit_all" -exec rm -f {} \;
endif

find . -name "core.*" -exec rm -f {} \; 
find . -name "fort.*" -exec rm -f {} \; 


rm running_linesearch
