/* Tell FFTW to behave very nicely and let us use normal operations for complex
 * arithmetic */
#include <complex.h>
#include <fftw3.h>
#include <math.h>

#include "config.h"
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
	complex double *in_rsp = arg_r->in_rsp;
	complex double *out_rsp = arg_r->out_rsp;

	/* numbers for second-order correction */

	complex double ***tidal_K = arg_r->tidal_K;
	complex double **lagrangian_s = arg_r->lagrangian_s;
	complex double **field_gradient = arg_r->field_gradient;

	/* position in input/output array */
	size_t l = arg_r->l;
	size_t m = arg_r->m;
	size_t n = arg_r->n;
	/* information about k-space */
	size_t X = arg_r->X;
	double real_spacing = arg_r->real_spacing;

	size_t idx = field_rsp_index(l, m, n, X);
	complex double field_val = in_rsp[idx];

	/* calculate correction */
	complex double delta_2 = 0.0;
	//eprintf("  %f + %fi   ", creal(field_val), cimag(field_val));
	
	delta_2 += 17.0 / 21.0 * field_val * field_val;
	for (size_t i = 0; i < 3; ++i) {
		for (size_t j = 0; j < 3; ++j) {
			delta_2 += 2.0 / 7.0 * tidal_K[i][j][idx] * tidal_K[i][j][idx];
		}
		complex double q =  lagrangian_s[i][idx] * field_gradient[i][idx];
		//eprintf("q = %f+%fi  ", creal(q), cimag(q));
		delta_2 -= q;
	}


	//eprintf("  %f + %fi   ", creal(delta_2), cimag(delta_2));
	out_rsp[idx] = field_val + delta_2;

	return 0;
}


void smooth(complex double *ksp, size_t KX, double mode_spacing)
{
	for (size_t l = 0; l < KX; ++l) {
		for (size_t m = 0; m < KX; ++m) {
			for (size_t n = 0; n < KX; ++n) {
				// Gaussian smoothing
				double k = index_to_k(l, m, n, KX, mode_spacing);
				ksp[field_rsp_index(l, m, n, KX)] *= (complex double) smoothing_gaussian(k);
			}
		}
	}
}

double smoothing_gaussian(double k)
{
	// in Mpc/h
	const double SMOOTH_LENGTH = PARAM_SMOOTH_LEN;
	const double SMOOTH_LENGTH_SQ = SMOOTH_LENGTH * SMOOTH_LENGTH;
	// follows from FT of a Gaussian in 3D
	// This Gaussian is normalised in real space, so not in reciprocal space.
	return exp(- k * k * SMOOTH_LENGTH_SQ / 2);
	//return (k > SMOOTH_LENGTH) ? 0.0 : 1.0;
}
