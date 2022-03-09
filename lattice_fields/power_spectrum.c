#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>

#include "util.h"
#include "power_spectrum.h"

double spec_flat(double k)
{
	/* just a uniform power spectrum for testing purposes */
	return 1.0;
}
double spec_gaussian(double k)
{
	return  exp(- k * k * pow(1.0e2, 2)) * 1e8;
}

double spec_delta(double k)
{
	if (k >= 1e-3)
		return 0.0;
	return 1.0e9;

}

double spec_linear(double k)
{
	return k;
}


static double bbks_f(double x);
static const double EPSILON = 1e-8;

double spec_bbks(double k)
{
	/* units of h/Mpc in k, (Mpc/h)^3 in A */
	const double A = 5.1e4;
	const double k_eq = 0.01;
	const double n = 0.967;

	return A * pow(k / k_eq, n) * pow(bbks_f(k / k_eq), 2);

}


static double bbks_f(double x)
{
	double Q = log(1 + 0.171 * x) / (0.171 * x + EPSILON);
	double R = 1 + 0.274 * x + pow(1.18 * x, 2) + pow(0.399 * x, 3) + pow(0.49 * x, 4);
	return Q * pow(R, -0.25);
}

/* TODO: write a more general field correlator function */
/* It will be a fairly simple copy-paste of this one */

void power_spectrum(complex double *field, size_t KX, double mode_spacing,
		double *k_buffer, double *bin_buffer, size_t *n_buffer, double k_min,
		double k_max, size_t n_bins)
{
	/* This is so we can find the mean k and power in each bin */
	//int *n_buffer = calloc(n_bins, sizeof(int));
	//actually, n_buffer is now provided to us
	
	/* For buffer assignment */
	double l_k_min = log(k_min);
	double l_k_max = log(k_max);
	double l_range = l_k_max - l_k_min;


	/* mode volume in k-space/density of states: if we get PS by doing
	 * an integral, equivalent in DFT is multiplying by V where it's nonzero*/
	double dV = pow(mode_spacing / (2 * M_PI), 3);

	/* Step through the field components in row-major order */
	for (size_t l = 0; l < KX; ++l) {
		for (size_t m = 0; m < KX; ++m){
			for (size_t n = 0; n < KX/2 + 1; ++n) {
				/* Find out what k is */
				double k = index_to_k(l, m, n, KX, mode_spacing);

				/* Find out which bin this k goes in. They are uniform
				 * in log-space. */
				double l_k = log(k);
				int bin = (int) trunc(n_bins * (l_k - l_k_min) / l_range);
				/* if k is out of the desired range, don't record it */

				if (bin < 0 || bin >= (int) n_bins)
				{
					continue;
				}

				/* Record the power, wavevector, and where it ended up */
				n_buffer[bin] += 1;

				complex double f_val = field[field_index(l, m, n, KX)];

				//complex double f_val_2 = field[field_index((KX - l) % KX, (KX - m) % KX, (KX -n) % KX, KX)];
				//complex double power = f_val * f_val_2;
				//eprintf("%f+%fi  ", creal(power), cimag(power));
				double power = f_val * conj(f_val) * dV;
				
				bin_buffer[bin] += power;

				k_buffer[bin] += k;
				/*
				if (k < 5.0) {
					eprintf("%f %f\n", k, f_imag);
				}
				*/

			}
		}
	}

	/* Step through the bins and take means of k and power */
	for (size_t i = 0; i < n_bins; ++i) {
		bin_buffer[i] /= n_buffer[i];
		k_buffer[i] /= n_buffer[i];
	}

	//free(n_buffer);
}




void correlator(complex double *field_1, complex double *field_2, size_t KX,
		double mode_spacing, double *k_buffer, double *bin_buffer,
		size_t *n_buffer, double k_min, double k_max, size_t n_bins)
{
	/* This is so we can find the mean k and power in each bin */
	//int *n_buffer = calloc(n_bins, sizeof(int));
	//actually, n_buffer is now provided to us
	
	/* For buffer assignment */
	double l_k_min = log(k_min);
	double l_k_max = log(k_max);
	double l_range = l_k_max - l_k_min;


	/* mode volume in k-space/density of states: if we get PS by doing
	 * an integral, equivalent in DFT is multiplying by V where it's nonzero*/
	double dV = pow(mode_spacing / (2 * M_PI), 3);

	/* Step through the field components in row-major order */
	for (size_t l = 0; l < KX; ++l) {
		for (size_t m = 0; m < KX; ++m){
			for (size_t n = 0; n < KX/2 + 1; ++n) {
				/* Find out what k is */
				double k = index_to_k(l, m, n, KX, mode_spacing);

				/* Find out which bin this k goes in. They are uniform
				 * in log-space. */
				double l_k = log(k);
				int bin = (int) trunc(n_bins * (l_k - l_k_min) / l_range);
				/* if k is out of the desired range, don't record it */

				if (bin < 0 || bin >= (int) n_bins)
				{
					continue;
				}

				/* Record the correlator, wavevector, and where it ended up */
				n_buffer[bin] += 1;

				complex double f1_val = field_1[field_index(l, m, n, KX)];
				complex double f2_val = field_2[field_index(l, m, n, KX)];

				double corr = f1_val * conj(f2_val) * dV;
				
				bin_buffer[bin] += corr;

				k_buffer[bin] += k;

			}
		}
	}

	/* Step through the bins and take means of k and correlator */
	for (size_t i = 0; i < n_bins; ++i) {
		bin_buffer[i] /= n_buffer[i];
		k_buffer[i] /= n_buffer[i];
	}

}
