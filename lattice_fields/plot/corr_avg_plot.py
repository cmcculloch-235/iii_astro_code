#!/usr/bin/env python3
import numpy as np
# Power spectrum

import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np


# load the data
assert(len(sys.argv) > 3)
power_data_file_base = sys.argv[1]
out_file_basename = sys.argv[2]
N_FILES = int(sys.argv[3])

# fourth argument specifies output format
EXTENSION = "pdf"
if len(sys.argv) > 4:
    EXTENSION = sys.argv[4]




# trick to initialise power_data
with open(power_data_file_base+ "0.dat") as power_data_fd:
    power_data_rows = np.array([[ float(j) for j in l.strip().split(" ")] for l in power_data_fd.readlines()[1:]])
    power_data = power_data_rows.transpose()

# iterate through data files and average power spectra.
for i in range(1, N_FILES):
    with open(power_data_file_base + str(i)+".dat") as power_data_fd:
        power_data_rows = np.array([[ float(j) for j in l.strip().split(" ")] for l in power_data_fd.readlines()[1:]])
        power_data += power_data_rows.transpose()

power_data /= N_FILES

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
stdev = power_data[4] * 1/np.sqrt(N_FILES * power_data[1])
axis.plot(power_data[2], power_data[4] + stdev, "--", color="red")
axis.plot(power_data[2], power_data[4] - stdev, "--", color="red")
axis.plot(power_data[2], power_data[4] + 2*stdev, "--", color="orange")
axis.plot(power_data[2], power_data[4] - 2*stdev, "--", color="orange")
#axis.plot(power_data[2], power_data[3], color="black", label="Linear")
axis.plot(power_data[2], power_data[3], color="green", label="$\\langle \\delta \\delta_\\tau \\rangle $")

axis.legend()

axis.set_xlabel("$k/ (h/\mathrm{Mpc})$")
axis.set_ylabel("$\\langle \\delta \\delta_\\tau \\rangle (k) / (\mathrm{Mpc}/h)^3$")
axis.set_title("Average of " + str(N_FILES) + " correlators")

axis.set_xscale("log", basex=10)
axis.set_yscale("log", basey=10)

axis.set_ylim(bottom=1e3)


plt.tight_layout()
plt.savefig(out_file_basename + "." + EXTENSION)

