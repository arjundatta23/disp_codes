#!/bin/bash

#python ~/PhD_Work/code/own/dispersion_codes/disp_fkMUSIC/music_supplement.py ${1} << END
#n
#y
#disp.bg.ray.gz
#0 3
#Aperture 2460 km
#END

inpdir=${1}
for file in ${inpdir}/*.pckl ; do
	python ~/PhD_Work/code/own/dispersion_codes/disp_fkMUSIC/music_supplement.py $file << END
n
n
END
	mv "music_pickle_plot.png" ${file/pckl/png}
done
