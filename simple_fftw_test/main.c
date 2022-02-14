#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <string.h>

#include "util.h"

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
	double L = 1000.0;
	double mode_spacing = 2.0 * M_PI / L;
	double real_spacing = L / X;

	/* Generate a Gaussian function */
	for (size_t i = 0; i < KX; ++i) {
		for (size_t j = 0; j < KX; ++j) {
			for (size_t k = 0; k < KX/2 + 1; ++k) {
				double q = index_to_k(i, 0, 0, KX, 1.0);
				double r = index_to_k(j, 0, 0, KX, 1.0);
				double s = index_to_k(k, 0, 0, KX, 1.0);

				double p = q * q + 2 * r*r + 3 * s*s ;

				size_t idx = field_index(i, j, k, KX);
				complex double f = cexp(-I * i);
				field_ksp_buf[idx] = f;
				size_t idx2 = field_index((KX-i)%KX, (KX-j)%KX, (KX-k)%KX, KX);
				field_ksp_buf[idx2] = conj(f);
				if (!i && !j && !k) {
					field_ksp_buf[idx2] = creal(f) * (1.0 );
				}
				if (i == j && j == k && k == KX/2) {
					field_ksp_buf[idx2] = creal(f);
				}

				//size_t idx2 = field_index((KX-i)%KX, (KX - j) % KX, (KX - k) % KX, KX);
				//field_ksp_buf[idx2] = exp( - q * q) * (1.0 - I);
				//if (j > 0) {
					//eprintf("%ld %ld %ld   " , l, m, n);
				//	field_ksp_buf[field_index(0, (KX - i) % KX, KX - j, KX)]
				//		= exp( - kx * kx - ky * ky);
				//}
			}
		}
	}
	//field_ksp_buf[KN - 1] = 100.0;
	//field_ksp_buf[1] = 100.0;

	for (size_t i = 0; i < X; ++i) {
		for (size_t j = 0; j < X; ++j) {
			//size_t idx = field_rsp_index(0, i, j, X);
			printf("%f+%fi  ", creal(field_ksp_buf[j + X * i]), cimag(field_ksp_buf[j + X * i]));
		}
		printf("\n");
	}
	fftw_execute(plan_k_to_r);

	printf("\n\n\n");
	//for (size_t i = 0; i < N; ++i) {
	//	printf("%ld %f+i%f\n", i, creal(field_rsp_buf[i]), cimag(field_rsp_buf[i]));
	//}
	for (size_t i = 0; i < X; ++i) {
		for (size_t j = 0; j < X; ++j) {
			size_t idx = field_rsp_index(i, j, 0, X);
			printf("%f+%fi  ", creal(field_rsp_buf[idx]), cimag(field_rsp_buf[idx]));
		}
		return 0;
		printf("\n");
	}



	/* ********************* *
	 * Clean up to be polite *
	 * ********************* */
	eprintf("Destroying plans...");
	fftw_destroy_plan(plan_r_to_k);
	fftw_destroy_plan(plan_k_to_r);

	eprintf("Freeing buffers...");
	fftw_free(field_ksp_buf);
	fftw_free(field_rsp_buf);
	eprintf("Done!\n");
	fftw_cleanup_threads();
	return 0;
}
