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
EXTENSION = "png"
if len(sys.argv) > 2:
	EXTENSION = sys.argv[2]

with open(data_file) as data_fd:
	data_rows = np.array([[ float(j) for j in l.strip().split(" ")] for l in data_fd.readlines()])
	data = data_rows.transpose()

# LaTeX mode
matplotlib.rcParams["text.usetex"] = True
matplotlib.rcParams["savefig.dpi"] = 300
matplotlib.rcParams["font.size"] = 18




# C211 correction
fig, axis = plt.subplots(1, 1)
axis.plot(data[0], data[1])

axis.set_xlabel("$k/ (h/\mathrm{Mpc})$")
axis.set_ylabel("$C^{211}_{(12)}(k) / (\mathrm{Mpc}/h)^3$")
axis.set_title("$C^{211}_{(12)}$ correction")

axis.set_xscale("log", basex=10)
axis.set_yscale("log", basey=10)

plt.tight_layout()
plt.savefig("out/delta_12/C211." + EXTENSION)


# C112 correction
fig, axis = plt.subplots(1, 1)
axis.plot(data[0], data[2])

axis.set_xlabel("$k/ (h/\mathrm{Mpc})$")
axis.set_ylabel("$C^{112}_{(12)}(k) / (\mathrm{Mpc}/h)^3$")
axis.set_title("$C^{112}_{(12)}$ correction")

axis.set_xscale("log", basex=10)
axis.set_yscale("log", basey=10)

plt.tight_layout()
plt.savefig("out/delta_12/C112." + EXTENSION)


# 1-loop correction
fig, axis = plt.subplots(1, 1)
axis.plot(data[0], data[3], color="blue")
axis.plot(data[0], -data[3], "--", color="blue")

axis.set_xlabel("$k/ (h/\mathrm{Mpc})$")
axis.set_ylabel("$C^{1L}_{(12)} / (\mathrm{Mpc}/h)^3$")
axis.set_title("Total 1-loop correlator")

axis.set_xscale("log", basex=10)
axis.set_yscale("log", basey=10)

plt.tight_layout()
plt.savefig("out/delta_12/1l_correction." + EXTENSION)



# Comparison of corrections
fig, axis = plt.subplots(1, 1)
axis.plot(data[0], data[1], color="blue", label="C211")
axis.plot(data[0], -data[1], "--", color="blue", label="")
axis.plot(data[0], data[2], color="orange", label="C112")
axis.plot(data[0], -data[2], "--", color="orange", label="")
axis.legend()

axis.set_xlabel("$k/ (h/\mathrm{Mpc})$")
axis.set_ylabel("$C_{(12)}(k) / (\mathrm{Mpc}/h)^3$")
axis.set_title("Individual 1-loop terms")

axis.set_xscale("log", basex=10)
axis.set_yscale("log", basey=10)

plt.tight_layout()
plt.savefig("out/delta_12/1l_correction_comparison." + EXTENSION)
