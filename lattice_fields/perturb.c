/* Tell FFTW to behave very nicely and let us use normal operations for complex
 * arithmetic */
#include <complex.h>
#include <fftw3.h>
#include <math.h>

#include "util.h"
#include "vector.h"

#include "perturb.h"

void *perturb_2(void *arg)
{
	/* right now strictly local but tidal effects will make this worth
	 * parallelising */
	/* Unpack the arguments */
	struct perturb_arg *arg_r = (struct perturb_arg *) arg;
	/* k-spaces of input and output fields */
	double *in_rsp = arg_r->in_rsp;
	double *out_rsp = arg_r->out_rsp;
	/* position in input/output array */
	size_t l = arg_r->l;
	size_t m = arg_r->m;
	size_t n = arg_r->n;
	/* information about k-space */
	size_t X = arg_r->X;
	double real_spacing = arg_r->real_spacing;

	/* calculate correction */
	double delta_2 = 0.0;

	delta_2 = 5.0 / 7.0 * pow(in_rsp[field_rsp_index(l, m, n, X)], 2);

	out_rsp[field_rsp_index(l, m, n, X)] += delta_2;

	return 0;
}


#ifdef AEIOU
/* this code could be useful but don't compile it */
static double F_2(double *k1, double *k2);
//const double EPSILON = 1e-4;
const double IR_CUTOFF = 1e-3;
static double F_2(double *k1, double *k2)
{
	/* Expects k1 and k2 to be vectors of length 3 */
	double k1k2 = vec_dot(k1, k2);
	double k1k1 = vec_dot(k1, k1);
	double k2k2 = vec_dot(k2, k2);
	if (k1k1 < pow(IR_CUTOFF, 2) || k2k2 < pow(IR_CUTOFF, 2)) {
		return 0.0;
	}
/*
	double A = 5.0/7.0;
	double B = 2.0 / 7.0 * pow(k1k2, 2) / (k1k1 * k2k2 + pow(EPSILON, 4));
	double C = k1k2 / 2 * (1 / (k1k1 + pow(EPSILON, 2)) + 1 / (k2k2 + pow(EPSILON, 2)));
*/

	double A = 5.0/7.0;
	double B = 2.0 / 7.0 * pow(k1k2, 2.0) / (k1k1 * k2k2);
	double C = k1k2 / 2 * (1 / (k1k1) + 1 / (k2k2));
	//printf("%f ", A+ B+ C);
	return A + B + C;
}







void *perturb_2(void *arg)
{
	/* Unpack the arguments */
	struct perturb_arg *arg_r = (struct perturb_arg *) arg;
	/* k-spaces of input and output fields */
	fftw_complex *in_ksp = arg_r->in_ksp;
	fftw_complex *out_ksp = arg_r->out_ksp;
	/* position in input/output array */
	size_t l = arg_r->l;
	size_t m = arg_r->m;
	size_t n = arg_r->n;
	/* information about k-space */
	size_t KX = arg_r->KX;
	double mode_spacing = arg_r->mode_spacing;


	/* calculate correction */
	fftw_complex delta_2 = 0.0;
	double dV = pow((double) KX, -3.0);

	double k[3];
	double k_1[3];
	double k_diff[3];

	index_to_vec_k(l, m, n, KX, mode_spacing, k);

	for (size_t p = 0; p < KX; ++p) {
		for (size_t q = 0; q < KX; ++q) {
			for (size_t r = 0; r < KX/2 + 1; ++r) {

				index_to_vec_k(p, q, r, KX, mode_spacing, k_1);
				//vec_sub(k, k_1, k_diff);

				/* Read off field values */
				complex double delta_1_1 = in_ksp[field_index(p, q, r, KX)];
				/* need to work out what idcs s, t, u correspond to k - k_1 */
				size_t s = k_diff_index(l, p, KX);
				size_t t = k_diff_index(m, q, KX);
				size_t u = k_diff_index(n, r, KX);
				index_to_vec_k(s, t, u, KX, mode_spacing, k_diff);
				complex double delta_1_2;

				/* if u is out of bounds for rfft, bring it back in
				 * using real FFT property */
				delta_1_2 = c2r_fft_access(s, t, u, KX, in_ksp);

				delta_2 += dV * F_2(k_1, k_diff) * delta_1_1 * delta_1_2;

				//eprintf("%f (%f %f %f) (%f %f %f) (%f %f %f) (%f %f %f) \n", KX * mode_spacing, k[0], k[1], k[2], k_1[0], k_1[1], k_1[2], k_diff[0], k_diff[1], k_diff[2], k[0] - k_1[0] - k_diff[0], k[1] - k_1[1] - k_diff[1], k[2] - k_1[2] - k_diff[2]);

				/* do modes at k_1 and -k_1 if r nonzero */
				if (r) {
					/* update the vectors but stay within first Brioullin zone */

					s = k_sum_index(l, p, KX);
                    t = k_sum_index(m, q, KX);
                    u = k_sum_index(n, r, KX);
					index_to_vec_k(s, t, u, KX, mode_spacing, k_diff);
					delta_1_2 = c2r_fft_access(s, t, u, KX, in_ksp);
					/* easier for first field because k_1 and -k_1 are
					 * certainly in 1BZ */
					vec_sca_mul(k_1, -1, k_1);
					/* update the field using real FFT property */
					delta_1_1 = conj(delta_1_1);
					/* more complex handling of the second field */
					/* find the indices for the field at k(l, m, n) +
					 * k_1(p, q, r) */

					delta_1_2 = c2r_fft_access(s, t, u, KX, in_ksp);

					delta_2 += dV * F_2(k_1, k_diff) * delta_1_1 * delta_1_2;
				}
			}
		}
	}

	/* write to output buffer */
	out_ksp[field_index(l, m, n, KX)] += delta_2;
	return 0;
}
#endif
