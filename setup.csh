#! /bin/csh

set directory=/home/jishnu/project/magnetic_inversions/working3/data 
set nmasterpixels=`wc -l < $directory/master.pixels`

## if (! -d $directory/tt/$itername) then
 ## mkdir $directory/tt/$itername
## endif

@ i = 1
while ($i <= $nmasterpixels) 
 
 if ($i < 10) then
  set forwardname=forward0$i
  set adjointname=adjoint0$i
 else
  set forwardname=forward$i
  set adjointname=adjoint$i
 endif

 if (! -d $directory/$forwardname) then
  mkdir $directory/$forwardname
 endif

 if (! -d $directory/$adjointname) then
  mkdir $directory/$adjointname
 endif
 @ i ++
end

 if (! -d $directory/update) then
  mkdir $directory/update
 endif

 if (! -d $directory/status) then
  mkdir $directory/status
 endif

 if (! -d $directory/kernel) then
  mkdir $directory/kernel
 endif

 if (! -d $directory/tt) then
  mkdir $directory/tt
 endif

 if (! -d $directory/data) then
  mkdir $directory/data
 endif



