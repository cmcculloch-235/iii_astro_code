#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <string.h>

#include "util.h"
#include "field_gen.h"
#include "power_spectrum.h"

int main(int argc, char *argv[])
{
	/* *********** *
	 * Set up FFTW * 
	 * *********** */
	/* Based on `info fftw` and Jeong's thesis */

	const int N_THREADS = 8;

	/* Documentation says to enable multi-threading before calling ANY fftw
	 * functions. Presumably inclides fftw_malloc. */
	fftw_init_threads();
	fftw_plan_with_nthreads(N_THREADS);
	

	/* Number of points per side of box in real space */
	size_t X = 256;

	/* Numbers of points in k-space */
	size_t KX = X;
	size_t KX_3 = X / 2 + 1;

	/* Number of modes, relevant for normalisation of FFTW output */
	size_t N = X * X * X;
	size_t KN = KX * KX * KX_3;

	/* k-space and real-space buffers */
	fftw_complex *field_ksp_buf = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * KN);
	double *field_rsp_buf = (double *) fftw_malloc(sizeof(double) * N);

	/* create plan for Fourier transforms */
	eprintf("Planning...");
	fftw_plan plan_r_to_k = fftw_plan_dft_r2c_3d(X, X, X, field_rsp_buf,
			field_ksp_buf, FFTW_MEASURE);


	/* create plan for inverse Fourier transforms */
	eprintf("Planning inverse...");
	fftw_plan plan_k_to_r = fftw_plan_dft_c2r_3d(KX, KX, KX, field_ksp_buf,
			field_rsp_buf, FFTW_MEASURE);
	eprintf("Done!\n");


	/* ************ *
	 * Set up Field * 
	 * ************ */

	/* Physical size of the box in units of Mpc/h */
	double L = 5000.0;
	double mode_spacing = 2.0 * M_PI / L;
	double real_spacing = L / X;
	
	/* Generate the linear field into the Fourier-space buffer */
	eprintf("Generating spectrum...");
	gen_field(field_ksp_buf, KX, mode_spacing, &spec_bbks);
	eprintf("Done!\n");


	/* ************ *
	 * Do some work *
	 * ************ */

	/* Inverse FFT the field */
	/* this overwrites the field buffer, so make a copy */
	fftw_complex *pre_fft = calloc(KN, sizeof(fftw_complex));
	memcpy(pre_fft, field_ksp_buf, KN * sizeof(fftw_complex));
	fftw_execute(plan_k_to_r);
	/* normalisation */
	for (size_t i = 0; i < N; ++i) {
		field_rsp_buf[i] /= N;
	}
	/* FFT the iFFTed field. Final field is in field_ksp_buf again*/
	fftw_execute(plan_r_to_k);



	/* *************** *
	 * Generate output *
	 * *************** */

	/* Extract the power spectrum before and after FFT and print it to 
	 * standard output */
	/*
	 * As a consistency check, make duplicate lists of ks and mode numbers in each bin
	 */
	size_t n_bins = 60;
	double *k_buffer_pre = calloc(n_bins, sizeof(double));
	size_t *n_buffer_pre = calloc(n_bins, sizeof(size_t));
	double *power_buffer_pre = calloc(n_bins, sizeof(double));
	double *k_buffer_post = calloc(n_bins, sizeof(double));
	size_t *n_buffer_post = calloc(n_bins, sizeof(size_t));
	double *power_buffer_post = calloc(n_bins, sizeof(double));

	power_spectrum(pre_fft, KX, mode_spacing, k_buffer_pre, power_buffer_pre,
			n_buffer_pre,
			// 5 mode_spaceing  up to 1.3KX mode_spacing / 2
			 5  * mode_spacing, 1.3 * KX / 2 * mode_spacing, n_bins);

	power_spectrum(field_ksp_buf, KX, mode_spacing, k_buffer_post,
			power_buffer_post, n_buffer_post, 5 * mode_spacing,
			1.3 * KX / 2 * mode_spacing, n_bins);

	/* print the power spectra to stdout */
	printf("bin n_pre n_post k_pre/h/Mpc k_post power_pre/(Mpc/h)^3 power_post power_reference\n");
	for (size_t i = 0; i < n_bins; ++i) {
		printf("%ld %ld %ld %f %f %f %f %f\n", i, n_buffer_pre[i], n_buffer_post[i],
				k_buffer_pre[i], k_buffer_post[i], power_buffer_pre[i], power_buffer_post[i],
				spec_bbks(k_buffer_pre[i]));
	}

	/* ******** *
	 * Clean up *
	 * ******** */
	eprintf("Destroying plans...");
	fftw_destroy_plan(plan_r_to_k);
	fftw_destroy_plan(plan_k_to_r);

	eprintf("Freeing buffers...");
	fftw_free(field_ksp_buf);
	fftw_free(field_rsp_buf);
	
	free(k_buffer_pre);
	free(k_buffer_post);
	free(n_buffer_pre);
	free(n_buffer_post);
	free(power_buffer_pre);
	free(power_buffer_post);

	fftw_free(pre_fft);
	eprintf("Done!\n");
	fftw_cleanup_threads();
	return 0;
}
