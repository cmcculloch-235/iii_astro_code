#!/bin/bash
./progress_plots.py ../data/1_loop.dat &&
	./progress_plots.py ../data/1_loop.dat eps &&
	cp out/progress/*.png ~/public_html/plots/progress/ &&
	cp out/progress/*.eps ~/public_html/plots/progress/eps/
