#!/bin/bash
./spectrum_plots.py ../data/1_loop.dat && ./spectrum_plots.py ../data/1_loop.dat eps && cp out/power_spectrum/*.png ~/public_html/plots/power_spectrum/ && cp out/power_spectrum/*.eps ~/public_html/plots/power_spectrum/eps/
