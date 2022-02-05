#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>


#include "util.h"


static double index_to_k_helper(size_t l, size_t KX, double mode_spacing);

int eprintf(const char *restrict format, ...)
{
	
	va_list ap;
	va_start(ap, format);
	int ret =  vfprintf(stderr, format, ap);
	va_end(ap);
	return ret;
}


double index_to_k(size_t l, size_t m, size_t n, size_t KX, double mode_spacing)
{
	double kl = index_to_k_helper(l, KX, mode_spacing);
	double km = index_to_k_helper(m, KX, mode_spacing);
	double kn = index_to_k_helper(n, KX, mode_spacing);


	return sqrt(kl * kl + km * km + kn * kn);
}

static double index_to_k_helper(size_t l, size_t KX, double mode_spacing)
{
	/* Accounts for the ordering of the modes in FFTW output */
	int mode_number = (int) l;
	if (l > KX / 2) 
		mode_number = ((int) l) - KX;
	
	return mode_spacing * mode_number;
}

void index_to_vec_k(size_t l, size_t m, size_t n, size_t KX,
		double mode_spacing, double *out)
{
	out[0] = index_to_k_helper(l, KX, mode_spacing);
	out[1] = index_to_k_helper(m, KX, mode_spacing);
	out[2] = index_to_k_helper(n, KX, mode_spacing);
}


static long signed int k_diff_helper(size_t l, size_t KX);

/* Given e.g. first indices of k and k_1, find first index of k - k_1 */
size_t k_diff_index(size_t l, size_t r, size_t KX)
{
	//eprintf("this function is broken");

	int q_1 = k_diff_helper(l, KX);
	int q_2 = k_diff_helper(r, KX);

	return (q_1 - q_2) % (signed long int) KX;

}

static long signed int k_diff_helper(size_t l, size_t KX)
{
	if (l <= KX / 2)
		return (long signed int) l;
	return (long signed int) l - KX;

}


size_t k_sum_index(size_t l, size_t r, size_t KX)
{
	int q_1 = k_diff_helper(l, KX);
	int q_2 = k_diff_helper(r, KX);

	return (q_1 + q_2) % (signed long int) KX;


}

fftw_complex c2r_fft_access(size_t l, size_t m, size_t n, size_t KX,
		fftw_complex *field)
{
	if (n <= KX/2) {
		return field[field_index(l, m, n, KX)];
	}
	/* n corresponded to a negative frequency. Use real FFT property */
	size_t p, q, r;
	p = k_diff_index(0, l, KX);
	q = k_diff_index(0, m, KX);
	r = k_diff_index(0, n, KX);
	return conj(field[field_index(p, q, r, KX)]);

}

/* Converts from [l][m][n] notation to [l + ()m + ()()n] */
size_t field_index(size_t l, size_t m, size_t n, size_t KX)
{
	return n + KX * m + KX * KX * l;
}
size_t field_rsp_index(size_t l, size_t m, size_t n, size_t X)
{
	return n + X * m + X * X * l;
}


