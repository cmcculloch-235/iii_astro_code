#!/bin/bash

./main > data/with_fft.dat
plot/plot.py data/with_fft.dat plot/with_fft
