#ifndef POWER_SPECTRUM_H_INC
#define POWER_SPECTRUM_H_INC

double spec_flat(double k);
double spec_linear(double k);
double spec_bbks(double k);

/* bin_buffer: mean power in the nth bin; k_buffer: the mean k for that bin */
void power_spectrum(fftw_complex *field, size_t KX, double mode_spacing, 
		double *k_buffer, double *bin_buffer, size_t *n_buffer, double k_min,
		double k_max, size_t n_bins);


#endif