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


void gen_tidal_K_ksp(void *in, void *out, size_t index, void *general_args);
void add_K_corr(void *in, void *out, size_t index, void *general_args);

struct gen_tidal_K_ksp_arg {
	size_t i;
	size_t j;
	double mode_spacing;
	double real_spacing;
	double KX;
};

double smoothing_gaussian(double k);
void smooth(complex double *ksp, size_t KX, double mode_spacing);

#endif
