#!/usr/bin/env python3

# TODO: get line type to change with change of sign

import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np


# load the data
assert(len(sys.argv) > 1)
data_file = sys.argv[1]

# second argument specifies output format
EXTENSION = "pdf"
if len(sys.argv) > 2:
	EXTENSION = sys.argv[2]

with open(data_file) as data_fd:
	data_rows = np.array([[ float(j) for j in l.strip().split(" ")] for l in data_fd.readlines()])
	data = data_rows.transpose()

# LaTeX mode
matplotlib.rcParams["text.usetex"] = True
matplotlib.rcParams["savefig.dpi"] = 300
matplotlib.rcParams["font.size"] = 18




# Ratio of correction to linear
fig, axis = plt.subplots(1, 1)
axis.plot(data[0], data[1] + data[2])
axis.plot(data[0], data[1])

axis.set_xlabel("$k/ (h/\mathrm{Mpc})$")
axis.set_ylabel("$P_{1L}/P_L / (\mathrm{Mpc}/h)^3$")
axis.set_title("Ratio of 1-loop correction to linear spectrum")

axis.set_xscale("log", basex=10)
axis.set_yscale("log", basey=10)
#axis.set_xlim(right=0.5)

plt.tight_layout()
plt.savefig("out/lattice_comparison/power_spectrum." + EXTENSION)

