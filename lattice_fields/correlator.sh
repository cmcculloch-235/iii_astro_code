#!/bin/bash

N_AVGS=10
RANGE=$(seq 0 $(($N_AVGS - 1)))
PARALLEL_PLOTS=7

# First, clear output dir
rm plot/out/corr_avg/*.pdf

for i in $RANGE; do
	echo $i:
	./main > data/corr_avg/$i.dat
done

echo -n "Plotting..."
parallel --jobs "$PARALLEL_PLOTS" plot/corr_plot.py "data/corr_avg/{}.dat" "plot/out/corr_avg/{}" ::: $RANGE

pdfunite plot/out/corr_avg/*.pdf plot/out/corr_avg/individual.pdf
echo "Done!"

# Now make the average plot
echo -n "Average plot..."
plot/corr_avg_plot.py "data/corr_avg/" plot/out/corr_avg/avg $N_AVGS
echo "Done!"
