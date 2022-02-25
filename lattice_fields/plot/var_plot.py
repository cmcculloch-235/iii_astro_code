#!/usr/bin/env python3
import numpy as np
import scipy.integrate as spi

import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np


# Auxiliary function for BBKS linear power spectrum
def bbks_f(x):
	EPSILON = 1e-6
	F_1 = np.log(1 + 0.171 * x) / (EPSILON + 0.171 * x)
	F_2 = 1 + 0.274 * x + (1.18 * x) ** 2 + (0.399 * x) ** 3 + (0.49 * x) ** 4
	return F_1 * np.power(F_2, -1/4)

# BBKS linear power spectrum
def spec_bbks(k):
	#print(k, end=" ")
	# k has units of h/Mpc as currently implemented
	A = 5.1e4 # (Mpc/h)^3
	#A = 1
	k_eq = 0.01 # h/Mpc
	n = 0.967
	return A * np.power(k / k_eq, n) * (bbks_f(k / k_eq) ** 2)

# Expected variance by doing the sum
N_POINTS = 256
L = 2500
k_0 = 2 * np.pi / L
k_min = (- (N_POINTS // 2) + 1) * k_0
k_max = (N_POINTS // 2 ) * k_0

kx, ky, kz = np.meshgrid(np.linspace(k_min, k_max, num = N_POINTS), np.linspace(k_min, k_max, num = N_POINTS), np.linspace(k_min, k_max, num = N_POINTS))
k_grid = np.sqrt(kx * kx + ky * ky + kz * kz)
discrete_spec = spec_bbks(k_grid)

def discrete_sigma(smooth_scale):
    smoothing = np.exp(-smooth_scale * smooth_scale * k_grid * k_grid)
    return np.sum(discrete_spec * smoothing * (k_0 / (2 * np.pi)) ** 3)

# Expected variance by numerical quadrature (on a sphere not a box)
def nquad_sigma(smooth_scale):
    integrand = lambda k: spec_bbks(k) * np.exp(-smooth_scale**2 * k**2) * k ** 2
    return (1 / (2 * np.pi ** 2)) * spi.nquad(integrand, [[0, k_max]])[0]



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
axis.plot(power_data[0], power_data[1], color="green", label="Lattice variance")
axis.plot(power_data[0], [discrete_sigma(s) for s in power_data[0]], "--", color="red", label="Exact expectation")
axis.plot(power_data[0], [nquad_sigma(s) for s in power_data[0]], "--", color="yellow", label="Numerical quadrature")

axis.legend()

axis.set_xlabel("Smoothing scale $\\Lambda$")
axis.set_ylabel("$\\sigma^2(\Lambda)$")
axis.set_title("Variance $\\sigma^2$ against smoothing scale $\\Lambda$")

#axis.set_ylim(bottom=1)

axis.set_xscale("log", basex=10)
axis.set_yscale("log", basey=10)

plt.tight_layout()
plt.savefig(out_file_basename + "." + EXTENSION)

