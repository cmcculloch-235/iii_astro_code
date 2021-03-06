#!/bin/bash

#SCALES="0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0 8.5 9.0 9.5 10.0"
SCALES="10.0 20.0 30.0 40.0 50.0 60.0 70.0 80.0 90.0 100.0 110.0 120.0 130.0 140.0 150.0 160.0"
PARALLEL_PLOTS=7

# First, clear output dir
rm plot/out/variance/linear/*.pdf
rm data/variance/linear/lin.dat

for S in $SCALES; do
	cat config.h | sed -e "s/#define PARAM_SMOOTH_LEN .*/#define PARAM_SMOOTH_LEN $S/g" > config_tmp.h
	mv config_tmp.h config.h
	make clean
	make -j8
	./main >> data/variance/linear/lin.dat
done



plot/var_plot.py data/variance/linear/lin.dat plot/out/variance/linear/lin
