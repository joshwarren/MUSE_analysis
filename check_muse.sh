#!/bin/bash

echo $1

if [ -z $2 ]
then
	echo 'Please enter the command as: check_muse.sh [galaxy] [opt_MCdir] [opptional dir]'
	exit	
fi

cd /Data/muse/analysis/$1/$2/MC/$3/


max=$( awk '{print $3}' /Data/muse/analysis/$1/$2/setup/voronoi_2d_binning_output.txt | sort -n | tail -1 )

i=0
while [[ $i -lt $max+1 ]]
do 
    ls ${i}.dat &> /dev/null || echo "     "$i / $max
    let i="i+1"
done
# done
