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
    power_data_rows = np.array([[ float(j) for j in l.strip().split(" ")] for l in power_data_fd.readlines()[1:]])
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
axis.plot(power_data[2], power_data[4], color="red", label="Linear power spectrum")
axis.plot(power_data[2], power_data[4] * (1 + 1/np.sqrt(power_data[1])), "--", color="red")
axis.plot(power_data[2], power_data[4] * (1 - 1/np.sqrt(power_data[1])), "--", color="red")
axis.plot(power_data[2], power_data[4] * (1 + 2/np.sqrt(power_data[1])), "--", color="orange")
axis.plot(power_data[2], power_data[4] * (1 - 2/np.sqrt(power_data[1])), "--", color="orange")
#axis.plot(power_data[2], power_data[3], color="black", label="Linear")
axis.plot(power_data[2], power_data[3], color="green", label="$\\langle \\delta \\delta_\\tau \\rangle $")
axis.legend()

axis.set_xlabel("$k/ (h/\mathrm{Mpc})$")
axis.set_ylabel("$\\langle \\delta \\delta_\\tau \\rangle (k) / (\mathrm{Mpc}/h)^3$")
#axis.set_title("Power spectrum")

#axis.set_ylim(bottom=1)

axis.set_xscale("log", basex=10)
axis.set_yscale("log", basey=10)

plt.tight_layout()
plt.savefig(out_file_basename + "." + EXTENSION)

