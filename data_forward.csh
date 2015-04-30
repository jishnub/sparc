#! /bin/csh
setenv MPI_TYPE_MAX 1280280

set directory = "/home/jishnu/project/magnetic_inversions/working3/data"
set nmasterpixels = `wc -l < "$directory/master.pixels"`

@ iter = 0
@ log = 1

while ($log == 1)
 
 @ iter++
 @ log = 0
 if (-e $directory/update/master.pixels.$iter) then
  @ log = 1 
 endif

end
 @ iter--

 if (! -d $directory/tt/data.$iter) then
  mkdir $directory/tt/data.$iter
 endif

find $directory -name "linesearch" -exec rm -f {} \; 
touch $directory/compute_data

@ i = 1
foreach i (`seq -f "%02g" 1 $nmasterpixels`)

  cp Spectral Instruction

  mpiexec -np 1 ./sparc > $directory/forward$i/data_forward_output
  find $directory/forward$i -name "*partial*" -exec rm -f {} \; 
  find $directory/forward$i -name "*full*" -exec rm -f {} \; 
  mv $directory/forward$i/vz_cc.fits $directory/forward$i/data.fits
  cp $directory/forward$i/data.fits $directory/tt/data.$iter/$i.fits
  cp $directory/forward$i/data.fits $directory/data/$i.fits

end

if ($i == $nmasterpixels) then
 find $directory/status -name "forward*" -exec rm -f {} \; 
 find $directory -name "compute_data" -exec rm -f {} \; 
 find $directory/kernel -name "misfit" -exec rm -f {} \; 
endif
