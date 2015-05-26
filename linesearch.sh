#!/bin/bash
#PBS -N  linesearch
#PBS -l nodes=2:ppn=24
#PBS -o  output-linesearch
#PBS -e  error-linesearch
#PBS -l walltime=12:00:00
cd $PBS_O_WORKDIR
echo $PBS_JOBID
export TERM=xterm


directory="/scratch/jishnu/magnetic/data"

find $directory -name "linesearch" -exec rm -f {} \; 
find $directory -name "compute_data" -exec rm -f {} \; 
find $directory -name "compute_synth" -exec rm -f {} \; 

iter=`find $directory/update -name 'linesearch_[0-9][0-9]'|wc -l`
itername=`printf "%02d" $iter`

touch $directory/linesearch
echo "Starting iterations at "`date`

for lin in `seq 1 5`
do
    linzpd=`printf "%02d" $lin`
    cp $directory/update/test_c_"$lin".fits  $directory/model_c_"$linzpd".fits
    cp $directory/update/test_vectorpsi_"$lin".fits  $directory/model_vectorpsi_"$linzpd".fits
done

/usr/local/bin/pbsdsh python $PBS_O_WORKDIR/linesearch.py

nmasterpixels=`wc -l < $directory/master.pixels`
for lin in `seq -f "%02g" 1 5`
do
    for src in `seq -f "%02g" 1 $((nmasterpixels))`
    do
        cat $directory/kernel/misfit_"$src"_"$lin" >> $directory/update/linesearch_$itername
        cat $directory/kernel/misfit_all_"$src"_"$lin" >> $directory/update/linesearch_all_$itername
        rm $directory/kernel/misfit_"$src"_"$lin"
        rm $directory/kernel/misfit_all_"$src"_"$lin"
    done
done


find $directory/update -name "tested*" -exec rm -f {} \; 
find $directory -name "update.fits" -exec rm -f {} \; 
find $directory -name "linesearch" -exec rm -f {} \; 
find $directory/status -name "forward*" -exec rm -f {} \;

find . -name "core.*" -exec rm -f {} \; 
find . -name "fort.*" -exec rm -f {} \; 

for lin in `seq 1 5`
do
    linzpd=`printf "%02d" $lin`
    rm $directory/model_c_"$linzpd".fits
    rm $directory/model_vectorpsi_"$linzpd".fits
done

echo "Finished at "`date`
