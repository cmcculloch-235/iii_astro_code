import numpy as np
import scipy.integrate as spi



# Naive convolution over k-space based on numerical quadrature
# Returns convolution, uncertainty in a NP vector
# Two functions f, g, and a kernel K. Returns value of the convolution at
# vector k.
# f, g : R^3 -> R
# K : R^3 x R^3 -> R
# Bounds of integration selected by...
# arbitrarily for now
# also this one is hilariously inefficient
def naive_stupid(f, g, K, k):
	vectorise = lambda q1, q2, q3: np.array([q1, q2, q3], np.float64)
	integrand = lambda q: K(q, k - q) * f(q) * g(k - q)
	integrand_wrap = lambda q1, q2, q3: integrand(vectorise(q1, q2, q3))
	bound = 3
	options={"limit": 80, "epsabs": -1, "epsrel": 1e-4}
	result = spi.nquad(integrand_wrap, [[-bound, bound], [-bound, bound],
		[-bound, bound]], opts=options)
	# it's okay to multiply the error as well, as scipy complains if round-off
	# error is important (which would invalidate the multiplication)
	# Actually, stopped returning the errors for now
	return 1 / (2 * np.pi) ** 3 * result[0]


# Still carrying some legacy quirkiness because I haven't devectorised the power
# spectrum function and maybe I could now. But this should be more efficient
# without changing the interface. So:
#
# TODO: write a new version with a more sensible interface (less function prep-
# aration required by caller)?
# Also consider squashing the range of k into a finite interval to remove
# questions of UV cutoff (in theory)
def naive_old_interface(f, g, K, k):
	# k = [k, 0, 0].
	q_range = [0, max([10, 1.5*k[0]])]
	q_range = [0, 10]
	mu_range = [-1, 1]
	options={"limit": 300, "epsabs": -1, "epsrel": 1e-11}
	
	# integrand set up to accept two vectors. Integrand_wrap converts the
	# coordinates on the domain of integration to one of these vectors.
	# the other is provided as k.
	to_vector = lambda q, mu: np.array([q * mu, q * np.sqrt(1 - mu**2), 0])
	integrand = lambda q: K(q, k - q) * f(q) * g(k - q)
	# incorporate reduced-dimensonality measure
	integrand_wrap = lambda q, mu: q**2 * integrand(to_vector(q, mu))

	# Pulled out a factor of 2pi by reducing dimensionality of integral
	result = 2 * np.pi * spi.nquad(integrand_wrap,
			[q_range, mu_range], opts=options)[0]
	# it's okay to multiply the error as well, as scipy complains if round-off
	# error is important (which would invalidate the multiplication)
	# Actually, stopped returning the errors for now
	return 1 / (2 * np.pi) ** 3 * result

# Isotropic integral over k-space of function f
def iso_kspace_integral(f, max_q):
	options={"limit": 200, "epsabs": -1, "epsrel": 1e-11}
	q_range = [0, max_q]
	integrand = lambda q: q**2 * f(q)
	result = 4 * np.pi * spi.nquad(integrand, [q_range], opts=options)[0]
	return (1/2 * np.pi)**3 * result

