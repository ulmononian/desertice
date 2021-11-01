#!/usr/bin/bash

# Run from iteration_archive dir

#export path='/global/cscratch1/sd/trhille/Thwaites_1km_GIA_ensemble/N3_01yr_PIGL/iteration_archive'

export files1=""
export files2=""
for d in `ls -1v`
do
echo $d
export files1="$files1 ${d}/iceload.nc "
export files2="$files2 ${d}/uplift_GIA.nc "

done

echo $files1
#ncrcat $files1 iceload_all.nc
# I had to manually run the command after pasting the echo output.  Not sure why.


#echo $files2
#ncrcat ${files2} uplift_all.nc


