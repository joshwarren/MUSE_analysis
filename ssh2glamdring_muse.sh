#!/bin/bash

# Bash script to ssh all necessary files to glamdring in one go.
reduced=/Data/muse/$1/$1.clipped.fits
analysis=/Data/muse/analysis/$1
analysis1=$analysis/templates.txt
analysis2=$analysis/voronoi_2d_binning_output_$2.txt
analysis3=$analysis/voronoi_2d_binning_output2_$2.txt
p_file=~/MUSE/analysis/params.txt
g_file=/Data/muse/analysis/galaxies.txt


scp $reduced warrenj@glamdring.physics.ox.ac.uk:muse_cubes/$1/
scp -r $analysis1 $analysis2 $analysis3 warrenj@glamdring.physics.ox.ac.uk:analysis_muse/$1/
scp $p_file warrenj@glamdring.physics.ox.ac.uk:
scp $g_file warrenj@glamdring.physics.ox.ac.uk:analysis_muse/