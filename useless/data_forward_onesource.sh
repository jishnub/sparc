#! /bin/bash
#export MPI_TYPE_MAX=1280280
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/users/apps/ciao-4.7/ciao-4.7/ots/lib/
directory="/scratch/jishnu/magnetic/data"
nmasterpixels=`wc -l < "$directory/master.pixels"`

#~ echo $directory
#~ echo $nmasterpixels
#~ echo $PBS_VNODENUM

if [[ ! -d $directory/tt/data ]]
then
    mkdir $directory/tt/data
fi

find $directory -name "linesearch" -exec rm -f {} \; 
touch $directory/compute_data

src=`printf "%02d" $((PBS_VNODENUM+1))`

cp Spectral Instruction$src

mpiexec -np 1 ./sparc > $directory/forward$src/data_forward_output
find $directory/forward$src -name "*partial*" -exec rm -f {} \; 
find $directory/forward$src -name "*full*" -exec rm -f {} \; 
mv $directory/forward$src/vz_cc.fits $directory/forward$src/data.fits
cp $directory/forward$src/data.fits $directory/tt/data/$src.fits
cp $directory/forward$src/data.fits $directory/data/$src.fits



if [[ "$PBS_VNODENUM" == $((nmasterpixels-1)) ]]
then
    find $directory/status -name "forward*" -exec rm -f {} \; 
    find $directory -name "compute_data" -exec rm -f {} \; 
    find $directory/kernel -name "misfit" -exec rm -f {} \; 
fi
