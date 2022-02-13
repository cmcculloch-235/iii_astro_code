#!/bin/bash
./delta_12_plots.py ../data/delta_12.dat && ./delta_12_plots.py ../data/delta_12.dat eps&& cp out/delta_12/*.png ~/public_html/plots/delta_12/ && cp out/delta_12/*.eps ~/public_html/plots/delta_12/eps/
