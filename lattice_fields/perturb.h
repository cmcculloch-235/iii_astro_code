#ifndef PERTURB_H_INC
#define PERTURB_H_INC

struct perturb_arg {
	/* real-spaces of input and output fields */
	double *in_rsp; double *out_rsp;
	/* position in input/output array */
	size_t l; size_t m; size_t n;
	/* information about k-space */
	size_t X; double real_spacing;
};

void *perturb_2(void *arg);

#endif