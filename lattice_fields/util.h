#ifndef UTIL_H_INC
#define UTIL_H_INC
#include <complex.h>
#include <fftw3.h>
int eprintf(const char *restrict format, ...);

double index_to_k(size_t l, size_t m, size_t n, size_t KX, double mode_spacing);
void index_to_vec_k(size_t l, size_t m, size_t n, size_t KX, double mode_spacing,
		double *out);

/* Given e.g. first indices of k and k_1, find first index of k - k_1 */
size_t k_diff_index(size_t l, size_t r, size_t KX);
size_t k_sum_index(size_t l, size_t r, size_t KX);
/* accepts 3 indices which would be good for a complex fft and an rfft field */
/* Only useful when doing mode coupling, which requires weird index manipulation*/
fftw_complex c2r_fft_access(size_t l, size_t m, size_t n, size_t KX,
		fftw_complex *field);

size_t field_index(size_t l, size_t m, size_t n, size_t KX);
size_t field_rsp_index(size_t l, size_t m, size_t n, size_t X);

#endif
