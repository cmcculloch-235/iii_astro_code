#!/bin/bash

N_SCALES=13
RANGE=$(seq 0 $(($N_SCALES - 1)))
SCALES=(002 004 006 008 010 020 040 060 080 100 120 140 160)
PARALLEL_PLOTS=7

# First, clear output dir
rm plot/out/smoothing/out/*.pdf

for i in $RANGE; do
	echo $i:
	S=${SCALES[$i]}
	cat config.h | sed -e "s/#define PARAM_SMOOTH_LEN .*/#define PARAM_SMOOTH_LEN $S.0/g" > config_tmp.h
	mv config_tmp.h config.h
	make clean
	make -j8
	./main > data/smoothing/$S.dat
done

echo -n "Plotting..."
parallel --jobs "$PARALLEL_PLOTS" plot/sample_plot.py "data/smoothing/{}.dat" "plot/out/smoothing/out/{}" ::: ${SCALES[*]}

pdfunite plot/out/smoothing/out/*.pdf plot/out/smoothing/out/individual.pdf
echo "Done!"

echo "Done!"
