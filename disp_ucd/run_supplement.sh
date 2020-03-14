#!/bin/bash

python ~/PhD_Work/code/own/dispersion_codes/disp_ucd/ucd_supplement.py ${1} << END
disp.bg.ray.gz
0 5
1
Aperture 2460 km
4
END
