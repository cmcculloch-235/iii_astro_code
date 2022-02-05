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


