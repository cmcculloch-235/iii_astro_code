import numpy as np

# Provide hard-coded perturbative kernels as found in appendix B.1 of
# Large-Scale Galaxy Bias, and in Cosmology with High-Redshift Galaxy Surveys

# F_n is for the density at nth perturbative order; G_n for velocity divergence

# F_n_m is a function used for calculating F_n.

# All kn arguments are vectors with three components.


# All top-level (i.e. F_n, G_n) kernels need to check that none of the
# denominators of any fractions are zero to avoid problems.
# Unfortunately, the number of vectors like k12, k123, k1...n, etc goes as 2^n
# Instead, just add a very tiny number to all the denominators to smooth them
# out
EPSILON = 1e-8
# Now with dimensionally consistent regularisation

def F_2(k1, k2):
	A = 2/7 * (np.dot(k1, k2)**2) / (EPSILON ** 2 + np.dot(k1, k1) * np.dot(k2, k2))
	B = (np.dot(k1, k2) / 2) * (1 / (np.dot(k1, k1) + EPSILON) +
			1 / (np.dot(k2, k2) + EPSILON))
	return 5/7 + A + B

def G_2(k1, k2):
	A = 4/7 * (np.dot(k1, k2) ** 2) / ( EPSILON**2 + np.dot(k1, k1) * np.dot(k2, k2))
	B = (np.dot(k1, k2) / 2) * (1 / (EPSILON + np.dot(k1, k1)) +
			1 / (EPSILON + np.dot(k2, k2)))
	return 3/7 + A + B



def F_3_1(k1, k2, k3):
	k23 = k2 + k3
	return np.dot(k1, k23) * G_2(k2, k3) / (EPSILON**2 + np.dot(k1, k1) * np.dot(k23, k23))

def F_3_2(k1, k2, k3):
	k12 = k1 + k2
	return k12 / (EPSILON + np.dot(k12, k12)) * G_2(k1, k2)

def F_3_3(k1, k2, k3):
	return k1 / (EPSILON + np.dot(k1, k1)) * F_2(k2, k3)

def F_3(k1, k2, k3):
	k123 = k1 + k2 + k3
	# there are more efficient ways of implementing the cyclic permutations,
	# but this is okay when there is only one third-order kernel.
	# possibly symmetrisation may be an obstacle to computing very high-order
	# kernels?
	A = 2 / 54 * np.dot(k123, k123) * ( F_3_1(k1, k2, k3) +
			F_3_1(k2, k3, k1) + F_3_1(k3, k1, k2))
	B = 7 / 54 * np.dot(k123, F_3_2(k1, k2, k3) + F_3_2(k2, k3, k1) +
		F_3_2(k3, k1, k2))
	C = 7 / 54 * np.dot(k123, F_3_3(k1, k2, k3) + F_3_3(k2, k3, k1) +
		F_3_3(k3, k1, k2))
	return A + B + C
