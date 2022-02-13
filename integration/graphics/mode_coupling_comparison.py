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





# Comparison of terms, and their sum
fig, axis = plt.subplots(1, 1)
axis.plot(data[0], data[1], color="blue", label="$P_L(k)$")
axis.plot(data[0], data[2], color="orange", label="$M(k)$")
axis.plot(data[0], -data[2], "--", color="orange", label="")
axis.plot(data[0], data[1] + data[2], color="green", label="$P_L(k) + M(k)$")
axis.legend()

axis.set_xlabel("$k/ (h/\\mathrm{Mpc})$")
axis.set_ylabel("Terms $/ (\\mathrm{Mpc}/h)^3$")
axis.set_title("Terms in Lyman-$\\alpha$ depth-\ndensity correlator")

axis.set_xscale("log", basex=10)
axis.set_yscale("log", basey=10)

plt.tight_layout()
plt.savefig("out/mode_coupling/mode_coupling_terms." + EXTENSION)


# Ratio of terms
fig, axis = plt.subplots(1, 1)
axis.plot(data[0], data[2] / data[1], color="blue", label="$M(k) / P_L(k)$")
axis.plot(data[0], -data[2] / data[1], "--", color="blue", label="")
axis.legend()

axis.set_xlabel("$k/ (h/\\mathrm{Mpc})$")
axis.set_ylabel("Ratio")
axis.set_title("Ratio of terms in Lyman-$\\alpha$ depth-\ndensity correlator")

axis.set_xscale("log", basex=10)
axis.set_yscale("log", basey=10)

plt.tight_layout()
plt.savefig("out/mode_coupling/mode_coupling_ratio." + EXTENSION)


