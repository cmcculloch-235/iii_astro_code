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
	complex double *in_rsp = arg_r->in_rsp;
	complex double *out_rsp = arg_r->out_rsp;
	/* position in input/output array */
	size_t l = arg_r->l;
	size_t m = arg_r->m;
	size_t n = arg_r->n;
	/* information about k-space */
	size_t X = arg_r->X;
	double real_spacing = arg_r->real_spacing;

	complex double field_val = in_rsp[field_rsp_index(l, m, n, X)];

	/* calculate correction */
	complex double delta_2 = 0.0;
	
	delta_2 = 5.0 / 7.0 * field_val * field_val;


	out_rsp[field_rsp_index(l, m, n, X)] += delta_2;

	return 0;
}

double smoothing_gaussian(double k);

void smooth(complex *ksp, size_t KX, double mode_spacing)
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
	const double h = 0.676;
	// in Mpc/h, so the number is the length in MPc
	const double SMOOTH_LENGTH = 50.0 / h;
	const double SMOOTH_LENGTH_SQ = SMOOTH_LENGTH * SMOOTH_LENGTH;
	// follows from FT of a Gaussian in 3D
	// This Gaussian is normalised in real space, so not in reciprocal space.
	return exp(- k * k * SMOOTH_LENGTH_SQ / 2);
}
