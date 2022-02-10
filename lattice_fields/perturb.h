#ifndef PERTURB_H_INC
#define PERTURB_H_INC

struct perturb_arg {
	/* real-spaces of input and output fields */
	complex double *in_rsp; complex double *out_rsp;
	/* quantities for second-order corrections */
	complex double ***tidal_K;
	complex double **lagrangian_s;
	complex double **field_gradient;
	/* position in input/output array */
	size_t l; size_t m; size_t n;
	/* information about k-space */
	size_t X; double real_spacing;
};

void *perturb_2(void *arg);

void smooth(complex *ksp, size_t KX, double mode_spacing);

#endif
