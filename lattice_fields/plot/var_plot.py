#!/usr/bin/env python3
import numpy as np
# Power spectrum

import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np


# load the data
assert(len(sys.argv) > 2)
#power_data_file = "../data/sample.dat"
power_data_file = sys.argv[1]
out_file_basename = sys.argv[2]

# second argument specifies output format
EXTENSION = "pdf"
if len(sys.argv) > 3:
    EXTENSION = sys.argv[3]

with open(power_data_file) as power_data_fd:
    power_data_rows = np.array([[ float(j) for j in l.strip().split(" ")] for l in power_data_fd.readlines()])
    power_data = power_data_rows.transpose()

# LaTeX mode
matplotlib.rcParams["text.usetex"] = True
matplotlib.rcParams["savefig.dpi"] = 300
matplotlib.rcParams["font.size"] = 12



# Comparison of spectra
#matplotlib.rcParams["font.size"] = 18
fig, axis = plt.subplots(1, 1)
#axis.set_aspect(0.3)
axis.margins(x=0, y=0.03)
axis.plot(power_data[0], power_data[1], color="green", label="Second-order")

#axis.legend()

axis.set_xlabel("Smoothing scale $\\Lambda$")
axis.set_ylabel("$\\sigma^2(\Lambda)$")
axis.set_title("Variance $\\sigma^2$ against smoothing scale $\\Lambda$")

#axis.set_ylim(bottom=1)

axis.set_xscale("log", basex=10)
axis.set_yscale("log", basey=10)

plt.tight_layout()
plt.savefig(out_file_basename + "." + EXTENSION)

