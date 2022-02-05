#include <stdio.h>
#include <stdarg.h>
#include <complex.h>
#include <math.h>
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




fftw_complex c2r_fft_access(size_t l, size_t m, size_t n, size_t KX,
		fftw_complex *field)
{
	if (n <= KX/2) {
		return field[field_index(l, m, n, KX)];
	}
	/* n corresponded to a negative frequency. */
	/* for now just return 0 to avoid a segfault */
	return 0.0;

}

/* Converts from [l][m][n] notation to [l + ()m + ()()n] */
size_t field_index(size_t l, size_t m, size_t n, size_t KX)
{
	return n + (KX ) * m + (KX) * KX * l;
}
size_t field_rsp_index(size_t l, size_t m, size_t n, size_t X)
{
	return n + X * m + X * X * l;
}


