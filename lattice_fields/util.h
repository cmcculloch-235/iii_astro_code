#ifndef UTIL_H_INC
#define UTIL_H_INC
#include <complex.h>
#include <fftw3.h>
int eprintf(const char *restrict format, ...);

double index_to_k(size_t l, size_t m, size_t n, size_t KX, double mode_spacing);
void index_to_vec_k(size_t l, size_t m, size_t n, size_t KX, double mode_spacing,
		double *out);

size_t field_index(size_t l, size_t m, size_t n, size_t KX);
size_t field_rsp_index(size_t l, size_t m, size_t n, size_t X);

complex double discrete_ksp_gradient(size_t l, size_t m, size_t n, size_t component,
		size_t KX, double real_spacing);
complex double discrete_ksp_laplacian(size_t l, size_t m, size_t n, size_t KX, double real_spacing);

#endif
