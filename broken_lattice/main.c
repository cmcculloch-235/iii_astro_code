#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <string.h>

#include "config.h"
#include "util.h"
#include "field_gen.h"
#include "power_spectrum.h"
#include "thread_pool.h"
#include "perturb.h"

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

	/* Number of modes, relevant for normalisation of FFTW output */
	size_t N = X * X * X;
	size_t KN = KX * KX * KX;

	/* k-space and real-space buffers */
	fftw_complex *field_ksp_buf = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * KN);
	fftw_complex *field_rsp_buf = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N);


	/* create plan for Fourier transforms */
	eprintf("Planning...");

	/* Check for FFTW wisdom */
	int wisdom_success = fftw_import_wisdom_from_filename("FFTW_wisdom");
	if (wisdom_success) {
		eprintf("Found FFTW wisdom...");
	}
	
	fftw_plan plan_r_to_k = fftw_plan_dft_3d(X, X, X, field_rsp_buf,
			field_ksp_buf, FFTW_FORWARD, FFTW_MEASURE);


	/* create plan for inverse Fourier transforms */
	eprintf("Planning inverse...");
	/* you would think the third argument would be KX_3 but it isn't */
	fftw_plan plan_k_to_r = fftw_plan_dft_3d(KX, KX, KX, field_ksp_buf,
			field_rsp_buf, FFTW_BACKWARD, FFTW_MEASURE);
	
	/* Always save the wisdom; it may be useful */
	fftw_export_wisdom_to_filename("FFTW_wisdom");
	eprintf("Saved FFTW wisdom...");

	eprintf("Done!\n");


	/* ************ *
	 * Set up Field * 
	 * ************ */

	/* Physical size of the box in units of Mpc/h */
	/* TODO: work out a sensible value for this */
	double L = 5000.0;
	double mode_spacing = 2.0 * M_PI / L;
	double real_spacing = L / X;
	
	/* Generate the linear field into the Fourier-space buffer */
	double (*spec_fn) (double) = &spec_bbks;
	eprintf("Generating spectrum...");
	gen_field(field_ksp_buf, KX, mode_spacing, spec_fn);
	eprintf("Done!\n");



	/* second-order corrections are calculated in real space */
	/* so inverse FFT the field */
	/* this overwrites the field buffer, so make a copy */
	fftw_complex *linear_ksp = calloc(KN, sizeof(fftw_complex));
	memcpy(linear_ksp, field_ksp_buf, KN * sizeof(fftw_complex));

	/* Add second-order correction */
	/* need smoothed real-space field */
	complex double *smoothed_rsp = calloc(N, sizeof(fftw_complex));
	eprintf("Prepare real-space...");
	eprintf("Smooth...");
	// Now smooth the kspace field before doing non-linear processing
	fftw_complex *smoothed_ksp = calloc(KN, sizeof(fftw_complex));
	memcpy(smoothed_ksp, field_ksp_buf, KN * sizeof(fftw_complex));
	smooth(smoothed_ksp, KX, mode_spacing);
/*
 * Checks that the field is the FFT of real data, which it is.
 */
	for (size_t l = 0; l < KX; ++ l) {
		for (size_t m = 0; m < KX; ++ m) {
			for (size_t n = 0; n < KX; ++ n) {
				size_t idx = field_index(l, m, n, KX);
				size_t l2 = (KX - l) % KX;
				size_t m2 = (KX - m) % KX;
				size_t n2 = (KX - n) % KX;
				size_t idx2 = field_index(l2, m2, n2, KX);
				complex double f = smoothed_ksp[idx] - smoothed_ksp[idx2];
				complex double g = smoothed_ksp[idx] + smoothed_ksp[idx2];
				if (cimag(g) != 0.0 ){
					eprintf("(%ld, %ld, %ld): %f+%fi  ",l, m, n, creal(f), cimag(g));
				}
			}
		}

	}
/**/

	memcpy(field_ksp_buf, smoothed_ksp, KN * sizeof(fftw_complex));

	eprintf("iFFT field...");
	fftw_execute(plan_k_to_r);

	eprintf("Normalise...");
	for (size_t i = 0; i < N; ++i) {
		field_rsp_buf[i] /= N;
	}
	memcpy(smoothed_rsp, field_rsp_buf, N * sizeof(complex double));


/*
 * Checks that the field is real
 */
	for (size_t l = 0; l < X; ++ l) {
		for (size_t m = 0; m < X; ++ m) {
			for (size_t n = 0; n < X; ++ n) {
				size_t idx = field_rsp_index(l, m, n, X);
				complex double f = smoothed_rsp[idx];
				const double tolerance = 5e-10;
				double ratio = fabs(cimag(f) / (creal(f) + tolerance));
				if (ratio > tolerance ){
					eprintf("(%ld, %ld, %ld): %f+%fi: %e  ",l, m, n, creal(f), cimag(f), ratio);
				}
			}
		}

	}
/**/

	/* Variance estimation */
#ifdef PARAM_VARIANCE
	/*in reciprocal space, it's just sum_k 1/KN PL(k) */
	double variance_ksp = 0.0;
	for (size_t i = 0; i < KN; ++i) {
		/* divide by KN = multiply by volume element in k-space */
		/* analogous to <delta(k) delta(k')>' = int_k' <delta(k) delta(k')> */
		variance_ksp += pow(cabs(smoothed_ksp[i]), 2) / KN;
	}
	/* and this factor of KN is from summing/integrating over k to find the
	 * variance from the power spectrum */
	variance_ksp /= KN;

	double variance_rsp = 0.0;
	for (size_t i = 0; i < N; ++i) {
		variance_rsp += pow(cabs(smoothed_rsp[i]), 2);
	}
	variance_rsp /= N;

	printf("%f %f %f\n", PARAM_SMOOTH_LEN, variance_ksp, variance_rsp);
	eprintf("%f %f %f\n", PARAM_SMOOTH_LEN, variance_ksp, variance_rsp);

#endif





	/* FFT the corrected spectrum into second_order_buffer*/
	//fftw_execute(plan_r_to_k);



	/* ********************* *
	 * Clean up to be polite *
	 * ********************* */
	eprintf("Destroying plans...");
	fftw_destroy_plan(plan_r_to_k);
	fftw_destroy_plan(plan_k_to_r);

	eprintf("Freeing buffers...");
	fftw_free(field_ksp_buf);
	fftw_free(field_rsp_buf);
	free(linear_ksp);
	free(smoothed_ksp);
	free(smoothed_rsp);

	eprintf("Done!\n");
	fftw_cleanup_threads();
	return 0;
}
