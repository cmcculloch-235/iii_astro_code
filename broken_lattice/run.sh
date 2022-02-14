#!/bin/bash

N_AVGS=1
RANGE=$(seq 0 $(($N_AVGS - 1)))
PARALLEL_PLOTS=7

# First, clear output dir
rm plot/out/avg/*.pdf

for i in $RANGE; do
	echo $i:
	./main > data/avg/$i.dat
done

echo -n "Plotting..."
#for i in $RANGE; do
	# Individual plots in parallel
#	plot/sample_plot.py data/avg/$i.dat plot/out/avg/$i &
#done
#wait
parallel --jobs "$PARALLEL_PLOTS" plot/sample_plot.py "data/avg/{}.dat" "plot/out/avg/{}" ::: $RANGE

pdfunite plot/out/avg/*.pdf plot/out/avg/individual.pdf
echo "Done!"

# Now make the average plot
echo -n "Average plot..."
plot/avg_plot.py "data/avg/" plot/out/avg/avg $N_AVGS
echo "Done!"
