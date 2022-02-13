#!/usr/bin/env python3
import numpy as np
# Power spectrum

import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np


# load the data
#assert(len(sys.argv) > 1)
power_data_file = "../data/1_loop.dat"

# second argument specifies output format
EXTENSION = "png"
if len(sys.argv) > 2:
	EXTENSION = sys.argv[2]

with open(power_data_file) as power_data_fd:
	power_data_rows = np.array([[ float(j) for j in l.strip().split(" ")] for l in power_data_fd.readlines()])
	power_data = power_data_rows.transpose()

# LaTeX mode
matplotlib.rcParams["text.usetex"] = True
matplotlib.rcParams["savefig.dpi"] = 300
matplotlib.rcParams["font.size"] = 12



# Comparison of spectra
matplotlib.rcParams["font.size"] = 18
fig, axis = plt.subplots(1, 1, figsize=(6, 6.5))
#axis.set_aspect(0.3)
axis.margins(x=0, y=0.03)
# individual terms
axis.plot(power_data[0], power_data[2], color="orange", label="$P^{22}$")
axis.plot(power_data[0], -power_data[2], "--", color="orange")
axis.plot(power_data[0], power_data[3], color="teal", label="$P^{13}$")
axis.plot(power_data[0], -power_data[3], "--", color="teal")
# corrected and uncorrected
axis.plot(power_data[0], power_data[1], color="green", label="Linear")
axis.plot(power_data[0], power_data[5], color="black", label="NLO")

axis.legend()

axis.set_xlabel("$k/ (h/\mathrm{Mpc})$")
axis.set_ylabel("$P(k) / (\mathrm{Mpc}/h)^3$")
axis.set_title("Contributions to the power spectrum")

axis.set_xscale("log", basex=10)
axis.set_yscale("log", basey=10)

plt.tight_layout()
plt.savefig("out/progress/power_spectrum_comparison." + EXTENSION)



# Composite operators on common axes
# C12 

# load the C12_data
#assert(len(sys.argv) > 1)
C12_data_file = "../data/delta_12.dat"

with open(C12_data_file) as C12_data_fd:
	C12_data_rows = np.array([[ float(j) for j in l.strip().split(" ")] for l in C12_data_fd.readlines()])
	C12_data = C12_data_rows.transpose()


# load the C22_data
C22_data_file = "../data/delta_22.dat"


with open(C22_data_file) as C22_data_fd:
	C22_data_rows = np.array([[ float(j) for j in l.strip().split(" ")] for l in C22_data_fd.readlines()])
	C22_data = C22_data_rows.transpose()

matplotlib.rcParams["font.size"] = 18
fig, axis = plt.subplots(1, 1, figsize=(6, 6.5))
#axis.set_aspect(0.5)
axis.margins(x=0, y=0.03)

axis.set_xscale("log", basex=10)
axis.set_yscale("log", basey=10)

axis.set_xlabel("$k/ (h/\mathrm{Mpc})$")

matplotlib.rcParams["font.size"] = 15
axis.set_ylabel("Correlators in $\langle\delta(\mathbf{k}) (\delta^2)(\mathbf{k'}) \\rangle$ \n and $\langle(\delta^2)(\mathbf{k}) (\delta^2)(\mathbf{k'}) \\rangle$ $/(\mathrm{Mpc}/h)^3$")

matplotlib.rcParams["font.size"] = 18
axis.set_title("One-loop terms in composite operators")


# Power spectrum-type plot: C112 with unit velocity dispersion
axis.plot(C12_data[0], C12_data[2], color="orange", label="$C_{(12)}^{112}$")
#axis.plot(C12_data[0], C12_data[2], color="orange")

# C211 correction
# the one of definite magnitude

axis.plot(C12_data[0], C12_data[1], color="blue", label="$C_{(12)}^{211}$")
axis.plot(C12_data[0], -C12_data[1], "--", color="blue")



# C22
# C1111 correction
#
axis.plot(C22_data[0], C22_data[1], color="green", label="$C_{(22)}^{1111}$")
axis.plot(C22_data[0], -C22_data[1], "--", color="green")

axis.margins(x=0, y=0.03)

#axis.set_yticks([10**9, 5 * 10**9])
#axis.set_yticklabels(["$10^9$", "$5\cdot10^{9}$"])

axis.legend()

plt.tight_layout()
plt.savefig("out/progress/composite." + EXTENSION)
