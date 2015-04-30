#PBS -S /bin/csh

# Job Name:
#PBS -N SYNTHETICS
##PBS -q devel

#PBS -l select=1:ncpus=1:model=ivy
#PBS -l walltime=04:00:00

 module add mpi-sgi/mpt.1.26
 setenv MPI_TYPE_MAX 1280280

 set directory = /nobackup/shanasog/classic
 set nmasterpixels = `wc -l < "$directory/master.pixels.synth"`

 ## COPYING SYNTH PIXELS TO MAIN PIXELS
 cp $directory/master.pixels $directory/master.pixels.temp
 cp $directory/master.pixels.synth $directory/master.pixels

 if (! -d $directory/syn) then
   mkdir $directory/syn
 endif

 ## COMPUTING DATA CORRESPONDING TO SYNTH PIXELS

 touch $directory/compute_data

 if (! -d $directory/syn/data) then
   mkdir $directory/syn/data
 endif

 @ i = 1
 while ($i <= $nmasterpixels) 

  cp Spectral Instruction

  if ($i < 10) then
   mpiexec -np 1 ./sparc > $directory/forward0$i/output
   find $directory/forward0$i -name "*partial*" -exec rm -f {} \; 
   find $directory/forward0$i -name "*full*" -exec rm -f {} \; 
   mv $directory/forward0$i/vz_cc.fits $directory/syn/data/0$i.fits
  else
   mpiexec -np 1 ./sparc > $directory/forward$i/output
   find $directory/forward$i -name "*partial*" -exec rm -f {} \; 
   find $directory/forward$i -name "*full*" -exec rm -f {} \; 
   mv $directory/forward$i/vz_cc.fits $directory/syn/data/$i.fits
  endif

  @ i++
 end
 find $directory -name "compute_data" -exec rm -f {} \; 

 ### NOW COMPUTING THE SYNTHETICS FOR ALL ITERATIONS
 ### FOR SYNTH PIXELS
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
 
 @ iter--

 echo $iter

 touch $directory/compute_synth
### touch $directory/linesearch
 @ j = 0
 while ($j <= $iter)

  if ($j < 10) then
   set index = 0$j
  else
   set index = $j
  endif

 if (! -e $directory/update/syn$index) then

##   echo $directory/update/model$index.fits $directory/model.fits
  find $directory/status -name "forward*" -exec rm -f {} \; 
  cp $directory/update/model_$index.fits $directory/model.fits

  set forwardname = $directory/syn/model$index
##  echo $forwardname
  if (! -d $forwardname) then
   mkdir $forwardname
  endif

  @ i = 1
  while ($i <= $nmasterpixels) 

   if ($i < 10) then
    set pix = 0$i
   else
    set pix = $i
   endif

   set forward = $directory/forward$pix

##   echo $forward
   cp Spectral Instruction
   mpiexec -np 1 ./sparc > $forward/output
   cp $forward/vz_cc.fits $forwardname/$pix.fits

   find $forward -name "*partial*" -exec rm -f {} \; 
   find $forward -name "*full*" -exec rm -f {} \; 

   @ i++

  end
  touch $directory/update/syn$index
##  echo "touching "$directory/update/tested$index
 endif

 @ j++

end

 find $directory -name "compute_synth" -exec rm -f {} \; 
## find $directory/status -name "forward*" -exec rm -f {} \;

 mv $directory/master.pixels.temp $directory/master.pixels

find . -name "core.*" -exec rm -f {} \; 
find . -name "fort.*" -exec rm -f {} \; 
