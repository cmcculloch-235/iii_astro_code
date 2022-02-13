import numpy as np

# Auxiliary function for BBKS linear power spectrum
def bbks_f(x):
	EPSILON = 1e-6
	F_1 = np.log(1 + 0.171 * x) / (EPSILON + 0.171 * x)
	F_2 = 1 + 0.274 * x + (1.18 * x) ** 2 + (0.399 * x) ** 3 + (0.49 * x) ** 4
	return F_1 * np.power(F_2, -1/4)

# BBKS linear power spectrum
def bbks(k):
	#print(k, end=" ")
	# k has units of h/Mpc as currently implemented
	A = 5.1e4 # (Mpc/h)^3
	#A = 1
	k_eq = 0.01 # h/Mpc
	n = 0.967
	return A * np.power(k / k_eq, n) * (bbks_f(k / k_eq) ** 2)
