#!/usr/bin/env python3
import numpy as np

import multiprocessing

from util import eprint

import spec_linear
import kernels
import convolve

eprint("Is this running nice-ly?")

# linear spectrum. Needs to accept vector arguments for ease of convolution
def PL_mm_vec(k):
	return spec_linear.bbks(np.sqrt(np.dot(k, k)))

# (22) contribution to mass-mass power spectrum
# use isotropy of the system to convert scalar k argument to vector in convolution
def P22_mm(k):
	return 2 * convolve.naive_old_interface(PL_mm_vec, PL_mm_vec,
		lambda p, q: kernels.F_2(p, q) ** 2, np.array([k, 0, 0], np.float64))

# looks a bit complicated just because I don't want to re-write convolve for this one
#P31_mm = lambda k: 3 * PL_mm_vec(k) * convolve.naive_old_interface(PL_mm_vec, lambda q: 1,
#		lambda p, q: kernels.F_3(p, -p, p-q), np.array([k, 0, 0], np.float64))

def P31_mm(k):
	return 3 * PL_mm_vec(k) * convolve.naive_old_interface(PL_mm_vec, lambda q: 1,
		lambda p, q: kernels.F_3(p, -p, p+q), np.array([k, 0, 0], np.float64))

# get ready for parallel processing 
# n.b. pool.map doesn't play nicely acting directly on lambdas
# that is, passing a lambda as an argument or having a lambda as the mapped
# function.
# problem is pickling the lambda


# basic pool worker, already obsolete
#def spec_1loop(k):
#	eprint(k, "P22")
#	p22_term = P22_mm(k)
#	eprint(k, "P31")
#	p31_term = P31_mm(k)
#	return (k, p22_term, p31_term)


# to make more efficient use of pooled resources
def pool_worker(arg):
	# Function is the function to be called (power spectrum term)
	# name is for debugging
	k, func, name = arg
	eprint(k, name)
	val = func(k)
	# it might not be necessary to return name, but it's useful to have it in
	# pool_worker for debugging
	return (k, val, name)


# usually -3, 0, 50 
ks = np.logspace(-3, 0, num=50)
args = [ (k, f[0], f[1]) for k in ks for f in [[PL_mm_vec, "PL"], [P22_mm, "P22"], [P31_mm, "P31"]]]

# mainly because I don't want to screw up and make infinite processes
if __name__ == "__main__":
	pool = multiprocessing.Pool(36)
	results = pool.map(pool_worker, args)

##################
# Prepare output #
###################
# Format: <k> <linear> <22> <31> <1l-correction> <overall spectrum>
# Don't output errors for now

# efficiently collect the calculations with the same k
N_terms = 3

# looks like [ [k, PL, P22, P31] ]
# but the way it's written is less clear
k_terms = [[results[i][0]] +  [results[i+j][1] for j in range(0, N_terms)]  
	for i in range(0, len(results), N_terms)]

# Prepare the output: sum some terms to get 1-loop correction to and corrected
# value of power spectrum

def correct_1_loop(contribs):
	# N.B. 1-loop correction is  P22 + 2P13.
	P22 = contribs[0]
	P13 = contribs[1]
	return P22 + 2 * P13
output = [k_term + [correct_1_loop(k_term[2:]), k_term[1] + correct_1_loop(k_term[2:])] for k_term in k_terms]

# print the output
for line in output:
	for entry in line:
		print(entry, end=" ")
	print("", end="\n")

