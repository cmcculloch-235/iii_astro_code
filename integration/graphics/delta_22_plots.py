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




# C1111 correction
fig, axis = plt.subplots(1, 1)
axis.plot(data[0], data[1])

axis.set_xlabel("$k/ (h/\mathrm{Mpc})$")
axis.set_ylabel("$C^{1111}_{(22)}(k) / (\mathrm{Mpc}/h)^3$")
axis.set_title("$C_{(22)}$, 1-loop")

axis.set_xscale("log", basex=10)
axis.set_yscale("log", basey=10)

plt.tight_layout()
plt.savefig("out/delta_22/C1111." + EXTENSION)




