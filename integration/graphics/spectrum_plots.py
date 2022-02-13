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


# Linear spectrum
fig, axis = plt.subplots(1, 1)
axis.plot(data[0], data[1])

axis.set_xlabel("$k/ (h/\mathrm{Mpc})$")
axis.set_ylabel("$PL(k) / (\mathrm{Mpc}/h)^3$")
axis.set_title("Linear power spectrum")

axis.set_xscale("log", basex=10)
axis.set_yscale("log", basey=10)

plt.tight_layout()
plt.savefig("out/power_spectrum/linear." + EXTENSION)


# P22 correction
fig, axis = plt.subplots(1, 1)
axis.plot(data[0], data[2])

axis.set_xlabel("$k/ (h/\mathrm{Mpc})$")
axis.set_ylabel("$P^{(22)}(k) / (\mathrm{Mpc}/h)^3$")
axis.set_title("$P^{(22)}$ correction")

axis.set_xscale("log", basex=10)
axis.set_yscale("log", basey=10)

plt.tight_layout()
plt.savefig("out/power_spectrum/P22." + EXTENSION)


# P31 correction
fig, axis = plt.subplots(1, 1)
axis.plot(data[0], -1 * data[3])

axis.set_xlabel("$k/ (h/\mathrm{Mpc})$")
axis.set_ylabel("$-P^{(31)}(k) / (\mathrm{Mpc}/h)^3$")
axis.set_title("$-1\\times P^{(31)}$ correction")

axis.set_xscale("log", basex=10)
axis.set_yscale("log", basey=10)

plt.tight_layout()
plt.savefig("out/power_spectrum/P31." + EXTENSION)


# 1-loop correction
fig, axis = plt.subplots(1, 1)
axis.plot(data[0], data[4], color="blue")
axis.plot(data[0], -data[4], "--", color="blue")

axis.set_xlabel("$k/ (h/\mathrm{Mpc})$")
axis.set_ylabel("$P_{1L} / (\mathrm{Mpc}/h)^3$")
axis.set_title("Total 1-loop correction")

axis.set_xscale("log", basex=10)
axis.set_yscale("log", basey=10)

plt.tight_layout()
plt.savefig("out/power_spectrum/1l_correction." + EXTENSION)


# Ratio of correction to linear
fig, axis = plt.subplots(1, 1)
axis.plot(data[0], data[4] / data[1])

axis.set_xlabel("$k/ (h/\mathrm{Mpc})$")
axis.set_ylabel("$P_{1L}/P_L / (\mathrm{Mpc}/h)^3$")
axis.set_title("Ratio of 1-loop correction to linear spectrum")

axis.set_xscale("log", basex=10)
axis.set_xlim(right=0.5)
axis.set_ylim(top=2, bottom=-0.1)

plt.tight_layout()
plt.savefig("out/power_spectrum/correction_ratio." + EXTENSION)

# Overall spectrum
fig, axis = plt.subplots(1, 1)
axis.plot(data[0], data[5])

axis.set_xlabel("$k/ (h/\mathrm{Mpc})$")
axis.set_ylabel("$P(k) / (\mathrm{Mpc}/h)^3$")
axis.set_title("1-loop corrected spectrum")

axis.set_xscale("log", basex=10)
axis.set_yscale("log", basey=10)

plt.tight_layout()
plt.savefig("out/power_spectrum/1l_spectrum." + EXTENSION)

# Comparison of spectra
fig, axis = plt.subplots(1, 1)
axis.plot(data[0], data[5], label="Corrected")
axis.plot(data[0], data[1], label="Linear")
axis.legend()

axis.set_xlabel("$k/ (h/\mathrm{Mpc})$")
axis.set_ylabel("$P(k) / (\mathrm{Mpc}/h)^3$")
axis.set_title("Linear and 1-loop corrected power spectra")

axis.set_xscale("log", basex=10)
axis.set_yscale("log", basey=10)

plt.tight_layout()
plt.savefig("out/power_spectrum/1l_comparison." + EXTENSION)

# Comparison of corrections
fig, axis = plt.subplots(1, 1)
axis.plot(data[0], np.abs(data[2]), label="P22")
axis.plot(data[0], np.abs(data[3]), label="-P13")
axis.legend()

axis.set_xlabel("$k/ (h/\mathrm{Mpc})$")
axis.set_ylabel("$P(k) / (\mathrm{Mpc}/h)^3$")
axis.set_title("Individual 1-loop terms")

axis.set_xscale("log", basex=10)
axis.set_yscale("log", basey=10)

plt.tight_layout()
plt.savefig("out/power_spectrum/1l_correction_comparison." + EXTENSION)
