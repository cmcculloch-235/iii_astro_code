#!/bin/bash

N_SCALES=5
RANGE=$(seq 0 $(($N_SCALES - 1)))
#SCALES=(0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0 8.5 9.0 9.5 10.0)
SCALES=(20 40 60 80 100)
PARALLEL_PLOTS=7

# First, clear output dir
#rm plot/out/avg/*.pdf

for i in $RANGE; do
	echo $i:
	S=${SCALES[$i]}
	cat config.h | sed -e "s/#define PARAM_SMOOTH_LEN .*/#define PARAM_SMOOTH_LEN $S/g" > config_tmp.h
	mv config_tmp.h config.h
	make clean
	make -j8
	./main > data/smoothing/$S.dat
done

echo -n "Plotting..."
parallel --jobs "$PARALLEL_PLOTS" plot/sample_plot.py "data/smoothing/{}.dat" "plot/out/smoothing/{}" ::: ${SCALES[*]}

pdfunite plot/out/smoothing/*.pdf plot/out/smoothing/individual.pdf
echo "Done!"

echo "Done!"
