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
	double (*spec_fn) (double) = &spec_bbks;
	eprintf("Generating spectrum...");
	gen_field(field_ksp_buf, KX, mode_spacing, spec_fn);
	eprintf("Done!\n");


// for debugging
/*
	size_t n_bins_debug = 60;
	double *k_buffer_debug = calloc(n_bins_debug, sizeof(double));
	double *power_buffer_debug = calloc(n_bins_debug, sizeof(double));
	size_t *n_buffer_debug = calloc(n_bins_debug, sizeof(size_t));

	double PSPEC_MIN_debug = 0.007;
	double PSPEC_MAX_debug = 0.2;
	power_spectrum(field_ksp_buf, KX, mode_spacing, k_buffer_debug, power_buffer_debug,
			n_buffer_debug,
			 PSPEC_MIN_debug, PSPEC_MAX_debug, n_bins_debug);

	printf("bin n k/h/Mpc lin_power/(Mpc/h)^3 2_power reference\n");
	for (size_t i = 0; i < n_bins_debug; ++i) {
		printf("%ld %ld %f %f %f\n", i, n_buffer_debug[i], k_buffer_debug[i], power_buffer_debug[i],
				spec_fn(k_buffer_debug[i]));
	}
*/


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
	memcpy(field_ksp_buf, smoothed_ksp, KN * sizeof(fftw_complex));

	eprintf("iFFT field...");
	fftw_execute(plan_k_to_r);
	eprintf("Normalise...");
	for (size_t i = 0; i < N; ++i) {
		field_rsp_buf[i] /= N;
	}
	memcpy(smoothed_rsp, field_rsp_buf, N * sizeof(complex double));

	/* ************************************ *
	 * Quantities for non-linear processing *
	 * ************************************ */
	/* Taking derivatives/dividing by k, so IR regularise */
	const double EPSILON=5e-3;
	/* Generate tidal tensor */
	eprintf("Tidal tensor...");
	/* Have to use nested pointers to make it easy to pass to perturb_2 */
	complex double ***tidal_K;
	
	
	tidal_K = calloc(3, sizeof(complex double **));
	for (size_t i = 0; i < 3; ++i) {
		tidal_K[i] = calloc(3, sizeof(complex double **));
		for (size_t j = 0; j <= i; ++j) {
			tidal_K[i][j] = calloc(N, sizeof(complex double));
		}
	}

	for (size_t i = 0; i < 3; ++i) {
		for (size_t j = 0; j <= i; ++j) {
			if (i != j) {
				tidal_K[j][i] = tidal_K[i][j];
			}
		}
	}

	/* Compute each index of the tidal tensor */
	eprintf("Tensor indices...");

	for (size_t l = 0; l < KX; ++ l) {
		for (size_t m = 0; m < KX; ++ m) {
			for (size_t n = 0; n < KX; ++ n) {
				size_t idx = field_index(l, m, n, KX);

				double vec_k[3];
				index_to_vec_k(l, m, n, KX, mode_spacing, vec_k);
				double sca_k = index_to_k(l, m, n, KX, mode_spacing);

				for (size_t i = 0; i < 3; ++i) {
					for (size_t j = 0; j <= i; ++j) {
						tidal_K[i][j][idx] = vec_k[i] * vec_k[j] * smoothed_ksp[idx] /
							(sca_k * sca_k + EPSILON * EPSILON);
						if (i == j) {
							tidal_K[i][j][idx] -= smoothed_ksp[idx] / 3;
						}
					}
				}
			}
		}
	}


	eprintf("iFFT tensor...");
	for (size_t i = 0; i < 3; ++i) {
		for (size_t j = 0; j <= i; ++j){
			memcpy(field_ksp_buf, tidal_K[i][j], KN * sizeof(complex double));
			fftw_execute(plan_k_to_r);
			memcpy(tidal_K[i][j], field_rsp_buf, N * sizeof(complex double));
		}
	}

	eprintf("and normalise...");
	for (size_t i = 0; i < 3; ++i) {
		for (size_t j = 0; j <= i; ++j){
			for (size_t l = 0; l < N; ++l) {
				tidal_K[i][j][l] /= N;
			}
		}
	}




	/* Real-space Lagrangian displaecment */
	eprintf("Lagrangian displacement...");
	complex double **lagrangian_s = calloc(3, sizeof(complex double *));
	for (size_t i = 0; i < 3; ++i) {
		lagrangian_s[i] = calloc(N, sizeof(complex double));
	}
	for (size_t l = 0; l < KX; ++ l) {
		for (size_t m = 0; m < KX; ++ m) {
			for (size_t n = 0; n < KX; ++ n) {
				double vec_k[3];
				index_to_vec_k(l, m, n, KX, mode_spacing, vec_k);
				double sca_k = index_to_k(l, m, n, KX, mode_spacing);
				size_t idx = field_index(l, m, n, KX);
				for (size_t i = 0; i < 3; ++i) {
					lagrangian_s[i][idx] = I * vec_k[i] * smoothed_ksp[idx] /
						(sca_k * sca_k + EPSILON * EPSILON);
				}
			}
		}
	}
	eprintf("iFFT vector...");
	for (size_t i = 0; i < 3; ++i) {
		memcpy(field_ksp_buf, lagrangian_s[i], KN * sizeof(complex double));
		fftw_execute(plan_k_to_r);
		memcpy(lagrangian_s[i], field_rsp_buf, N * sizeof(complex double));
	}

	eprintf("and normalise...");
	for (size_t l = 0; l < N; ++l) {
		for (size_t i = 0; i < 3; ++i) {
			lagrangian_s[i][l] /= N;
		}
	}
	/* Real-space field gradient) */
	eprintf("Field gradient...");
	complex double **field_gradient = calloc(3, sizeof(complex double *));
	for (size_t i = 0; i < 3; ++i) {
		field_gradient[i] = calloc(N, sizeof(complex double));
	}
	for (size_t l = 0; l < KX; ++ l) {
		for (size_t m = 0; m < KX; ++ m) {
			for (size_t n = 0; n < KX; ++ n) {
				double vec_k[3];
				index_to_vec_k(l, m, n, KX, mode_spacing, vec_k);
				size_t idx = field_index(l, m, n, KX);
				for (size_t i = 0; i < 3; ++i) {
					field_gradient[i][idx] = I * vec_k[i] * smoothed_ksp[idx];
				}
			}
		}
	}
	eprintf("iFFT vector...");
	for (size_t i = 0; i < 3; ++i) {
		memcpy(field_ksp_buf, field_gradient[i], KN * sizeof(complex double));
		fftw_execute(plan_k_to_r);
		memcpy(field_gradient[i], field_rsp_buf, N * sizeof(complex double));
	}

	eprintf("and normalise...");
	for (size_t l = 0; l < N; ++l) {
		for (size_t i = 0; i < 3; ++i) {
			field_gradient[i][l] /= N;
		}
	}

	eprintf("Done!\n");

	


	struct perturb_arg *arg_buffer = calloc(N, sizeof(struct perturb_arg));
	struct pool_job *job_buffer = calloc(N, sizeof(struct pool_job));
	eprintf("Preparing second-order jobs...");
	for (size_t l = 0; l < X; ++l){
		for (size_t m = 0; m < X; ++m) {
			for (size_t n = 0; n < X; ++n) {
				size_t i = field_rsp_index(l, m, n, X);

				arg_buffer[i].in_rsp = smoothed_rsp;
				arg_buffer[i].out_rsp = field_rsp_buf;

				arg_buffer[i].tidal_K = tidal_K;
				arg_buffer[i].lagrangian_s = lagrangian_s;
				arg_buffer[i].field_gradient = field_gradient;


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
	struct pool_work work = {.jobs = job_buffer, .n_jobs = N, .n_collect = 4096};
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
	free(linear_ksp);
	free(smoothed_ksp);
	free(smoothed_rsp);

	for (size_t i = 0; i < 3; ++i) {
		for (size_t j = 0; j <= i; ++j) {
			free(tidal_K[i][j]);
		}
		free(lagrangian_s[i]);
		free(field_gradient[i]);
		free(tidal_K[i]);
	}
	free(lagrangian_s);
	free(field_gradient);
	free(tidal_K);
	eprintf("Done!\n");
	fftw_cleanup_threads();
	return 0;
}
