#!/usr/bin/env python3
import numpy as np
from numpy.random import default_rng
# Power spectrum

import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np


assert(len(sys.argv) > 1)
out_file_basename = sys.argv[1]

# second argument specifies output format
EXTENSION = "pdf"
if len(sys.argv) > 2:
    EXTENSION = sys.argv[2]


# Generate data for Taylor expanded and exact functions of a Gaussian
# random variable.
# Use the same realisation of the variable for each calculation
# Array content: (sigma, mean_exact, mean_taylor_1, mean_taylor_3, mean_taylor_5)

alpha = 1.7
def f_exact(d):
    return np.power(1 + d, alpha)

def falling_factorial(a, n):
    ret = 1
    for i in range(0, n):
        ret *= a - i
    return ret

def factorial(n):
    return falling_factorial(n, n)

def f_taylor_1(d):
    return 1 + d * alpha

def f_taylor_2(d):
    return 1 + d * alpha + d ** 2 * alpha * (alpha - 1) / 2
    
def f_taylor_3(d):
    return (f_taylor_1(d) + d**2 * falling_factorial(alpha, 2) / factorial(2) +
            d**3 * falling_factorial(alpha, 3) / factorial(6))

def f_taylor_5(d):
    return (f_taylor_3(d) + d**4 * falling_factorial(alpha, 4) / factorial(4) +
            d**5 * falling_factorial(alpha, 5) / factorial(5))

def f_taylor_10(d):
    return (f_taylor_5(d) +
            d**6 * falling_factorial(alpha, 6) / factorial(6) +
            d**7 * falling_factorial(alpha, 7) / factorial(7) +
            d**8 * falling_factorial(alpha, 8) / factorial(8) +
            d**9 * falling_factorial(alpha, 9) / factorial(9) +
            d**10 * falling_factorial(alpha, 10) / factorial(10))

def f_taylor_k(d, k):
    # kth order Taylor expansion
    ret = 0
    for i in range(0, k + 1):
        ret += d ** i * falling_factorial(alpha, i) / factorial(i)
    return ret


data = []
sigmas = np.geomspace(0.1, 10, num=100)
N_samples = 5000

print("Generating...")
# set up rng
rng = default_rng()
for s in sigmas:
    print(s)
    exact = 0
    taylor_orders = [1, 2, 3, 5, 30]
    n_orders = len(taylor_orders)
    taylors = [0 for i in range(0, n_orders)]
    for i in range(0, N_samples):
        d = -2
        while(d < - 1):
            d = rng.normal(scale=s)
        exact += f_exact(d)
        for j in range(0, n_orders):
            taylors[j] += f_taylor_k(d, taylor_orders[j])
    exact /= N_samples
    taylors = [t / N_samples for t in taylors]
    data.append([s, exact] + taylors)
# arrange the data into columns
data = np.array(data)
data = data.transpose()
print("Done!")
#print(data)




# LaTeX mode
matplotlib.rcParams["text.usetex"] = True
matplotlib.rcParams["savefig.dpi"] = 300
matplotlib.rcParams["font.size"] = 12


# Data format:
# bin n_pre n_post k_pre/h/Mpc k_post power_pre/(Mpc/h)^3 power_post power_reference

# Comparison of spectra
#matplotlib.rcParams["font.size"] = 18
fig, axis = plt.subplots(1, 1)
#axis.set_aspect(0.3)
axis.margins(x=0, y=0.03)
axis.plot(data[0], data[6], color="black", label="Thirtieth-order")
axis.plot(data[0], data[5], color="red", label="Quintic")
axis.plot(data[0], data[4], color="purple", label="Cubic")
axis.plot(data[0], data[3], linestyle="dotted", color="violet", label="Quadratic")
axis.plot(data[0], data[2], color="blue", label="Linear")
axis.plot(data[0], data[1], color="green", label="Exact")

axis.legend()

axis.set_xlabel("Standard deviation $\\sigma$ of $\\delta$")
axis.set_ylabel("$\\langle(1 + \\delta)^{1.7}\\rangle$")
axis.set_title("Monte-Carlo estimation of $\\langle(1 + \delta)^{1.7}\\rangle$ using\nexact and Taylor-exapanded functions")

axis.set_xscale("log", basex=10)
axis.set_yscale("log", basey=10)
axis.set_ylim(0.1, 100)

plt.tight_layout()
plt.savefig(out_file_basename + "." + EXTENSION)


