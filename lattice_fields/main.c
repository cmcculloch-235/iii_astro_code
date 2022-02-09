#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <string.h>

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
	double (*spec_fn) (double) = &spec_flat;
	eprintf("Generating spectrum...");
	gen_field(field_ksp_buf, KX, mode_spacing, spec_fn);
	eprintf("Done!\n");



	/* Add second-order correction */
	/* need somewhere to put them */
	complex double *quad_rsp_buffer = calloc(N, sizeof(fftw_complex));
	/* second-order corrections are calculated in real space */
	/* so inverse FFT the field */
	/* this overwrites the field buffer, so make a copy */
	fftw_complex *linear_ksp = calloc(KN, sizeof(fftw_complex));
	memcpy(linear_ksp, field_ksp_buf, KN * sizeof(fftw_complex));

	eprintf("Smooth...");
	// Now smooth the kspace field before doing non-linear processing
	smooth(field_ksp_buf, KX, mode_spacing);
	eprintf("R2K FFT...");
	fftw_execute(plan_k_to_r);
	eprintf("Normalise...");
	for (size_t i = 0; i < N; ++i) {
		field_rsp_buf[i] /= N;
	}
	eprintf("Done!");
	//fftw_execute(plan_r_to_k);

	
	memcpy(quad_rsp_buffer, field_rsp_buf, N * sizeof(complex double));


	struct perturb_arg *arg_buffer = calloc(N, sizeof(struct perturb_arg));
	struct pool_job *job_buffer = calloc(N, sizeof(struct pool_job));
	eprintf("Preparing second-order jobs...");
	for (size_t l = 0; l < X; ++l){
		for (size_t m = 0; m < X; ++m) {
			for (size_t n = 0; n < X; ++n) {
				size_t i = field_rsp_index(l, m, n, X);

				arg_buffer[i].in_rsp = quad_rsp_buffer;
				arg_buffer[i].out_rsp = field_rsp_buf;

				arg_buffer[i].l = l;
				arg_buffer[i].m = m;
				arg_buffer[i].n = n;
				arg_buffer[i].X = X;
				arg_buffer[i].real_spacing = real_spacing;

				job_buffer[i].function = &perturb_2;
				job_buffer[i].arg = (void *) &(arg_buffer[i]);
			}
		}
	}
	eprintf("Done!\n");
	struct pool_work work = {.jobs = job_buffer, .n_jobs = N, .n_collect = 1024};
	pool_run(&work,N_THREADS);

	free(arg_buffer);
	free(job_buffer);


	/* FFT the corrected spectrum into second_order_buffer*/
	fftw_execute(plan_r_to_k);



	/* Extract the power spectrum at first and second order and print it to 
	 * standard output */
	size_t n_bins = 60;
	double *k_buffer = calloc(n_bins, sizeof(double));
	double *power_buffer_lin = calloc(n_bins, sizeof(double));
	size_t *n_buffer = calloc(n_bins, sizeof(size_t));
	double *power_buffer_2 = calloc(n_bins, sizeof(double));

	double PSPEC_MIN = 0.007;
	double PSPEC_MAX = 0.2;

	power_spectrum(linear_ksp, KX, mode_spacing, k_buffer, power_buffer_lin,
			n_buffer,
			// 5 mode_spaceing  up to 1.3KX mode_spacing / 2
			 PSPEC_MIN, PSPEC_MAX, n_bins);

	/* zero out the k and n buffers for re-use */
	bzero(k_buffer, n_bins * sizeof(double));
	bzero(n_buffer, n_bins * sizeof(size_t));
	
	power_spectrum(field_ksp_buf, KX, mode_spacing, k_buffer,
			power_buffer_2, n_buffer, PSPEC_MIN, PSPEC_MAX, n_bins);

	/* print the power spectra to stdout */
	printf("bin n k/h/Mpc lin_power/(Mpc/h)^3 2_power reference\n");
	for (size_t i = 0; i < n_bins; ++i) {
		printf("%ld %ld %f %f %f %f\n", i, n_buffer[i], k_buffer[i], power_buffer_lin[i], power_buffer_2[i],
				spec_fn(k_buffer[i]));
	}

	/* ********************* *
	 * Clean up to be polite *
	 * ********************* */
	free(k_buffer);
	free(power_buffer_lin);
	free(n_buffer);
	free(power_buffer_2);
	eprintf("Destroying plans...");
	fftw_destroy_plan(plan_r_to_k);
	fftw_destroy_plan(plan_k_to_r);

	eprintf("Freeing buffers...");
	fftw_free(field_ksp_buf);
	fftw_free(field_rsp_buf);
	fftw_free(linear_ksp);
	free(quad_rsp_buffer);
	eprintf("Done!\n");
	fftw_cleanup_threads();
	return 0;
}
