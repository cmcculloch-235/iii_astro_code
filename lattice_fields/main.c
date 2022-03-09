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
#include "config.h"

#include "transformations.h"

int main(int argc, char *argv[])
{




	/* *********** *
	 * Set up FFTW * 
	 * *********** */
	/* Based on `info fftw` and Jeong's thesis */

	const int N_THREADS = 8;

	/*
	// thread map test
	size_t TEST_N = 300;
	struct test_arg arg = {.n = TEST_N};
	int *numbers = calloc(TEST_N, sizeof(int));
	int *new_numbers = calloc(TEST_N, sizeof(int));
	thread_map(map_test_function, &arg, numbers, sizeof(int), new_numbers,
			sizeof(int), TEST_N, N_THREADS);
	for (int i = 0; i < TEST_N; ++i){
		eprintf("%d ", new_numbers[i]);
	}
	return 0;
	*/




	/* Documentation says to enable multi-threading before calling ANY fftw
	 * functions. Presumably inclides fftw_malloc. */
	fftw_init_threads();
	fftw_plan_with_nthreads(N_THREADS);
	

	/* Number of points per side of box in real space */
	size_t X = 128;

	/* Numbers of points in k-space */
	size_t KX = X;

	/* Number of modes, relevant for normalisation of FFTW output */
	size_t N = X * X * X;
	size_t KN = KX * KX * KX;

	/* k-space and real-space buffers */
	complex double *field_ksp_buf = (complex double *) fftw_malloc(sizeof(complex double) * KN);
	complex double *field_rsp_buf = (complex double *) fftw_malloc(sizeof(complex double) * N);


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
	double L = 2500.0;
	double mode_spacing = 2.0 * M_PI / L;
	double real_spacing = L / X;
	// 1/(this * N) is also (deltak/2pi)^3
	double real_dV = real_spacing * real_spacing * real_spacing;
	
	/* Generate the linear field into the Fourier-space buffer */
	double (*spec_fn) (double) = &spec_bbks;
	eprintf("Generating spectrum...");
	/* might want to parallelise this, but managing the RNG is a little subtle */
	gen_field(field_ksp_buf, KX, mode_spacing, spec_fn);


	eprintf("Done!\n");

	fftw_execute(plan_k_to_r);

	/* FFT and inverse for debugging purposes: verify that they compose to 1 */
	struct normalisation_arg nn_arg = {.real_dV = real_dV, .N = N};
	thread_map(&normalise_k_to_r, (void *) &nn_arg,
			(void *) field_rsp_buf, sizeof(complex double),
			(void *) field_rsp_buf, sizeof(complex double),
			N, N_THREADS);
	fftw_execute(plan_r_to_k);
	thread_map(&normalise_r_to_k, (void *) &nn_arg,
			(void *) field_ksp_buf, sizeof(complex double),
			(void *) field_ksp_buf, sizeof(complex double),
			N, N_THREADS);



	/* second-order corrections are calculated in real space */
	/* so inverse FFT the field */
	/* this overwrites the field buffer, so make a copy */
	complex double *linear_ksp = calloc(KN, sizeof(complex double));
	memcpy(linear_ksp, field_ksp_buf, KN * sizeof(complex double));

	/* Add second-order correction */
	/* need smoothed real-space field */
	complex double *smoothed_rsp = calloc(N, sizeof(complex double));
	eprintf("Prepare real-space...");
	// Now smooth the kspace field before doing non-linear processing
	complex double *smoothed_ksp = calloc(KN, sizeof(complex double));
	memcpy(smoothed_ksp, field_ksp_buf, KN * sizeof(complex double));
#if (PARAM_SMOOTH)
	eprintf("Smooth...");
	smooth(smoothed_ksp, KX, mode_spacing);
#endif


	memcpy(field_ksp_buf, smoothed_ksp, KN * sizeof(complex double));

	eprintf("iFFT field...");
	fftw_execute(plan_k_to_r);

	eprintf("Normalise...");
	/*
	for (size_t i = 0; i < N; ++i) {
		field_rsp_buf[i] /= N * real_dV;
	}
	*/


	struct normalisation_arg n_arg = {.real_dV = real_dV, .N = N};
	thread_map(&normalise_k_to_r, (void *) &n_arg,
			(void *) field_rsp_buf, sizeof(complex double),
			(void *) field_rsp_buf, sizeof(complex double),
			N, N_THREADS);

	memcpy(smoothed_rsp, field_rsp_buf, N * sizeof(complex double));



	/* Variance estimation */
#if (PARAM_VARIANCE)
	eprintf("variance...");
	/*in reciprocal space, it's just sum_k 1/KN PL(k) */
	double variance_ksp = 0.0;
	for (size_t i = 0; i < KN; ++i) {
		/* divide by KN = multiply by volume element in k-space */
		/* analogous to <delta(k) delta(k')>' = int_k' <delta(k) delta(k')> */
		variance_ksp += pow(cabs(smoothed_ksp[i]), 2) / (KN * real_dV);
	}
	/* and this factor of KN is from summing/integrating over k to find the
	 * variance from the power spectrum */
	variance_ksp /= (KN * real_dV);

	double variance_rsp = 0.0;
	for (size_t i = 0; i < N; ++i) {
		variance_rsp += pow(cabs(smoothed_rsp[i]), 2);
	}
	variance_rsp /= N;

	printf("%f %f %f\n", PARAM_SMOOTH_LEN, variance_ksp, variance_rsp);
	eprintf("%f %f %f\n", PARAM_SMOOTH_LEN, variance_ksp, variance_rsp);
	return 0;
#endif







	/* ************************************ *
	 * Quantities for non-linear processing *
	 * ************************************ */

	complex double *nl_rsp_correction = calloc(N, sizeof(complex double));
	/* Laplacian is small at low k, so IR regularise */
	const double EPSILON=1e-6;
	/* Generate tidal tensor */
	eprintf("Tidal tensor...");

	for (size_t i = 0; i < 3; ++i) {
		for (size_t j = 0; j <= i; ++j) {
			/* Calculate the correction due to this component of the
			 * tidal tensor and immediately add it to nl_corrections
			 * Factors of 2 are stored in add_K_corr to avoid running expensive
			 * computations twice
			 */
			struct gen_tidal_K_ksp_arg gen_K_arg;

			gen_K_arg.i = i;
			gen_K_arg.j = j;
			gen_K_arg.mode_spacing = mode_spacing;
			gen_K_arg.real_spacing = real_spacing;
			gen_K_arg.KX = KX;

			thread_map(&gen_tidal_K_ksp, (void *) &gen_K_arg, (void *) smoothed_ksp,
					sizeof(complex double), (void *) field_ksp_buf, sizeof(complex double),
					KN, N_THREADS);

			fftw_execute(plan_k_to_r);

			/* n_arg is defined above */
			thread_map(&normalise_k_to_r, (void *) &n_arg,
					(void *) field_rsp_buf, sizeof(complex double),
					(void *) field_rsp_buf, sizeof(complex double),
					N, N_THREADS);


			thread_map(&add_K_corr, (void *) &gen_K_arg, (void *) field_rsp_buf,
					sizeof(complex double), (void *) nl_rsp_correction, sizeof(complex double),
					N, N_THREADS);

		}
	}

	eprintf("Quadratic term...");
	thread_map(&add_quad_corr, NULL, (void *) smoothed_rsp,
			sizeof(complex double), (void *) nl_rsp_correction, sizeof(complex double),
			N, N_THREADS);

	/* This one is a little trickier to calculate using map, as it needs s.grad d,
	 * but just carry grad d in general_args and look it up using index */


	eprintf("Displacement-gradient term...");
	for (size_t i = 0; i < 3; ++i) {
		complex double *field_gradient = calloc(N, sizeof(complex double));
		complex double *field_displacement = calloc(N, sizeof(complex double));

		struct gen_ksp_grad_dis_arg gen_arg;

		gen_arg.mode_spacing = mode_spacing;
		gen_arg.real_spacing = real_spacing;
		gen_arg.KX = KX;
		gen_arg.i = i;

		/* Get the field gradient and send it to real space */
		thread_map(&gen_ksp_gradient, (void *) &gen_arg, (void *) smoothed_ksp,
				sizeof(complex double), (void *) field_ksp_buf, sizeof(complex double),
				KN, N_THREADS);
		fftw_execute(plan_k_to_r);
		thread_map(&normalise_k_to_r, (void *) &n_arg,
				(void *) field_rsp_buf, sizeof(complex double),
				(void *) field_rsp_buf, sizeof(complex double),
				N, N_THREADS);
		memcpy(field_gradient, field_rsp_buf, N * sizeof(complex double));

		/* Get the field displacement and send it to real space */
		thread_map(&gen_ksp_displacement, (void *) &gen_arg, (void *) smoothed_ksp,
				sizeof(complex double), (void *) field_ksp_buf, sizeof(complex double),
				KN, N_THREADS);
		fftw_execute(plan_k_to_r);
		thread_map(&normalise_k_to_r, (void *) &n_arg,
				(void *) field_rsp_buf, sizeof(complex double),
				(void *) field_rsp_buf, sizeof(complex double),
				N, N_THREADS);
		memcpy(field_displacement, field_rsp_buf, N * sizeof(complex double));

		/* Multiply the two and subtract from correction */
		struct add_ksp_grad_dis_arg add_arg;

		add_arg.field_gradient = field_gradient;

		thread_map(&add_grad_dis_corr, (void *) &add_arg,
				(void *) field_displacement, sizeof(complex double),
				(void *) nl_rsp_correction, sizeof(complex double),
				N, N_THREADS);

		free(field_gradient);
		free(field_displacement);

	}


	/* move correction into fftw buffer */
	memcpy(field_rsp_buf, nl_rsp_correction, N * sizeof(complex double));

	/* FFT and normalise */
	eprintf("FFT correction...");
	fftw_execute(plan_r_to_k);

	thread_map(&normalise_r_to_k, (void *) &n_arg,
			(void *) field_ksp_buf, sizeof(complex double),
			(void *) field_ksp_buf, sizeof(complex double),
			KN, N_THREADS);


	/* add it to the initial field */
	eprintf("Add to original...");
	
	thread_map(&add_nl_correction, NULL,
			(void *) smoothed_ksp, sizeof(complex double),
			(void *) field_ksp_buf, sizeof(complex double),
			KN, N_THREADS);

	complex double *nl_ksp = calloc(KN, sizeof(complex double));
	memcpy(nl_ksp, field_ksp_buf, KN * sizeof(complex double));

	//free(nl_rsp_correction);
	eprintf("Done!");


	size_t n_bins = 60;
	double K_MIN = 0.007;
	double K_MAX = 0.2;

#if PARAM_PSPEC_MODE
	/* Extract the power spectrum at first and second order and print it to 
	 * standard output */
	double *k_buffer = calloc(n_bins, sizeof(double));
	double *power_buffer_lin = calloc(n_bins, sizeof(double));
	size_t *n_buffer = calloc(n_bins, sizeof(size_t));
	double *power_buffer_2 = calloc(n_bins, sizeof(double));


//	power_spectrum(linear_ksp, KX, mode_spacing, k_buffer, power_buffer_lin,
	power_spectrum(smoothed_ksp, KX, mode_spacing, k_buffer, power_buffer_lin,
			n_buffer,
			// 5 mode_spaceing  up to 1.3KX mode_spacing / 2
			 K_MIN, K_MAX, n_bins);

	/* zero out the k and n buffers for re-use */
	bzero(k_buffer, n_bins * sizeof(double));
	bzero(n_buffer, n_bins * sizeof(size_t));
	
	power_spectrum(nl_ksp, KX, mode_spacing, k_buffer,
			power_buffer_2, n_buffer, K_MIN, K_MAX, n_bins);



	/* print the power spectra to stdout */
	printf("bin n k/h/Mpc lin_power/(Mpc/h)^3 2_power reference\n");
	for (size_t i = 0; i < n_bins; ++i) {
		printf("%ld %ld %f %f %f %f\n", i, n_buffer[i], k_buffer[i], power_buffer_lin[i], power_buffer_2[i],
				pow(smoothing_gaussian(k_buffer[i]), 2) * spec_fn(k_buffer[i]));
				//spec_fn(k_buffer[i]));
	}
	free(k_buffer);
	free(power_buffer_lin);
	free(n_buffer);
	free(power_buffer_2);

	/* free these here because they are freed earlier in the else below */
	free(linear_ksp);
	free(smoothed_ksp);
	free(smoothed_rsp);

#else

	/* Free these early because we're about to allocate more buffers */
	free(linear_ksp);
	free(smoothed_ksp);

	/* Generate the two fields to be correlated, using the non-linear field */
	/* as it stands, nl correction to rsp representation is in field_rsp_buf;
	 * therefore, add the original smoothed field to it*/

	thread_map(&add_nl_correction, NULL,
			(void *) smoothed_rsp, sizeof(complex double),
			(void *) field_rsp_buf, sizeof(complex double),
			N, N_THREADS);
	/* don't need smoothed_rsp anymore and we're about to allocate another
	 * buffer */
	free(smoothed_rsp);

	complex double *field_k_1 = calloc(N, sizeof(complex double));
	complex double *field_k_2 = calloc(N, sizeof(complex double));

	eprintf("Field transformations...");
	/* transformations are in real space */

	/* There are enough buffers to avoid allocating another for a copy of the
	 * real field */
	memcpy(field_k_1, field_rsp_buf, N * sizeof(complex double));
	
	/*
	thread_map(&PARAM_CORR_F1, NULL,
			(void *) field_rsp_buf, sizeof(complex double),
			(void *) field_rsp_buf, sizeof(complex double),
			KN, N_THREADS);
	*/
	PARAM_CORR_F1(field_rsp_buf, N, N_THREADS);

	/* FFT and normalise */
	fftw_execute(plan_r_to_k);
	thread_map(&normalise_r_to_k, (void *) &n_arg,
			(void *) field_ksp_buf, sizeof(complex double),
			(void *) field_ksp_buf, sizeof(complex double),
			KN, N_THREADS);

	/* restore copy of non-linear field */
	memcpy(field_rsp_buf, field_k_1, N * sizeof(complex double));
	/* copy transformed field to field_k_1 buffer */
	memcpy(field_k_1, field_ksp_buf, KN * sizeof(complex double));


	/* Do the same thing for the second field. No need to copy the non-linear
	 * field this time */
	PARAM_CORR_F2(field_rsp_buf, N, N_THREADS);

	fftw_execute(plan_r_to_k);
	thread_map(&normalise_r_to_k, (void *) &n_arg,
			(void *) field_ksp_buf, sizeof(complex double),
			(void *) field_ksp_buf, sizeof(complex double),
			KN, N_THREADS);

	/* copy transformed field to field_k_1 buffer */
	memcpy(field_k_2, field_ksp_buf, KN * sizeof(complex double));

	eprintf("Extract correlator...");

	/* Now extract the correlator */
	double *k_buffer = calloc(n_bins, sizeof(double));
	size_t *n_buffer = calloc(n_bins, sizeof(size_t));
	double *corr_buffer = calloc(n_bins, sizeof(double));

	correlator(field_k_1, field_k_1, KX, mode_spacing, k_buffer, corr_buffer,
			n_buffer, K_MIN, K_MAX, n_bins);

	/* print the correlator to stdout */
	printf("bin n k/h/Mpc corr/(Mpc/h)^3 PL\n");
	for (size_t i = 0; i < n_bins; ++i) {
		printf("%ld %ld %f %f %f\n", i, n_buffer[i], k_buffer[i], corr_buffer[i],
				pow(smoothing_gaussian(k_buffer[i]), 2) * spec_fn(k_buffer[i]));
	}
	eprintf("Done!\n");
	free(k_buffer);
	free(n_buffer);
	free(corr_buffer);
	free(field_k_1);
	free(field_k_2);

#endif

	/* ********************* *
	 * Clean up to be polite *
	 * ********************* */
	eprintf("Destroying plans...");
	fftw_destroy_plan(plan_r_to_k);
	fftw_destroy_plan(plan_k_to_r);

	eprintf("Freeing buffers...");
	fftw_free(field_ksp_buf);
	fftw_free(field_rsp_buf);
	free(nl_ksp);

	eprintf("Done!\n");
	fftw_cleanup_threads();
	return 0;
}
