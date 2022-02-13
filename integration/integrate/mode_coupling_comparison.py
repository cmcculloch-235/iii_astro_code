#!/usr/bin/env python3
import numpy as np

import multiprocessing

from util import eprint

import spec_linear
import kernels
import convolve

N_THREADS = 8

eprint("Is this running nice-ly?")

# linear spectrum. Needs to accept vector arguments for ease of convolution
def PL_mm_vec(k):
	return spec_linear.bbks(np.sqrt(np.dot(k, k)))


# Mode-coupling term
def M(k):
	return convolve.naive_old_interface(PL_mm_vec, PL_mm_vec, lambda p, q: kernels.F_2(p, q), np.array([k, 0, 0], np.float64))



# get ready for parallel processing 
# n.b. pool.map doesn't play nicely acting directly on lambdas
# that is, passing a lambda as an argument or having a lambda as the mapped
# function.
# problem is pickling the lambda


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
ks = np.logspace(-3, 0, num=100)
args = [ (k, f[0], f[1]) for k in ks for f in [[PL_mm_vec, "P_L"], [M, "M"]]]

# mainly because I don't want to screw up and make infinite processes
if __name__ == "__main__":
	pool = multiprocessing.Pool(N_THREADS)
	results = pool.map(pool_worker, args)

##################
# Prepare output #
###################
# Format: <k>  <P_L> <M> 
# Don't output errors for now

# efficiently collect the calculations with the same k
N_terms = 2

# looks like [ [k, P_L, M] ]
# but the way it's written is less clear
k_terms = [[results[i][0]] +  [results[i+j][1] for j in range(0, N_terms)]  
	for i in range(0, len(results), N_terms)]

# Prepare the output: sum some terms to get 1-loop correction to and corrected
# value of power spectrum

output = k_terms

# print the output
for line in output:
	for entry in line:
		print(entry, end=" ")
	print("", end="\n")

