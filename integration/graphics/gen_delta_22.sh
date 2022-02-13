#!/bin/bash
./delta_22_plots.py ../data/delta_22.dat && ./delta_22_plots.py ../data/delta_22.dat eps&& cp out/delta_22/*.png ~/public_html/plots/delta_22/ && cp out/delta_22/*.eps ~/public_html/plots/delta_22/eps/
