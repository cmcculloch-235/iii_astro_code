#include <math.h>
#include <complex.h>
#include <fftw3.h>
/* gsl has a good interface for Gaussian random numbers */
#include <gsl/gsl_randist.h>
/* to seed the RNG */
#include <time.h>

#include "util.h"
#include "field_gen.h"


void gen_field(fftw_complex *field_buffer, size_t KX, double mode_spacing,
		double (*power_spectrum) (double))
{
	/* Volume element in k-space is just 1/number of points in real space */
	double dV = pow((double) KX, -3);

	/* Set up the RNG. Type is luxury random numbers, which have good
	 * decorrelation.*/
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_ranlxd2);
	/* seed the RNG. Wants a 32-bit float*/
	gsl_rng_set(rng, time(NULL));
	
	/* Step through the modes, get the power spectrum, and generate the field */

	for (size_t l = 0; l < KX; ++l) {
		for (size_t m = 0; m < KX; ++m) {
			for (size_t n = 0; n < KX/2 + 1; ++n) {

				double k = index_to_k(l, m, n, KX, mode_spacing);
				/* evaluate power spectrum at k, then draw real and imaginary
				 * parts from corresponding Gaussian distributions */
				double ps = power_spectrum(k);
				double sigma = sqrt(ps / (2.0 * dV));
				//sigma = sqrt(ps /2.0);

				double real_part = gsl_ran_gaussian(rng, sigma);
				double imag_part = gsl_ran_gaussian(rng, sigma);



				/* we are stepping through in row-major order */
				field_buffer[field_index(l, m, n, KX)] = real_part + I * imag_part;
				
				// enforce real-ness of field
				field_buffer[field_index((KX - l) % KX, (KX - m) % KX, (KX - n) % KX, KX)]
					= real_part - I * imag_part;
				// special cases where -k and k modes are equal
				if ((l == 0 || (l == KX/2 && ! (KX % 2))) &&
						(m == 0 || (m == KX/2 && ! (KX % 2))) &&
						(n == 0 || (n == KX/2 && ! (KX % 2))) ){
					field_buffer[field_index(l, m, n, KX)] = real_part;
				}



				/* ********* *
				 * DEBUGGING *
				 * ********* **/
				#ifdef FG_DEBUG
				if (!(l ||m ||n)) {
					eprintf("Debugging field mode");
				}
				/*
				if (( l == KX/4 ) ) {
					real_part = KX * KX * KX * KX;
				}*/
				complex double field_val = cexp(-k * k - KX/2 * I * l) * KX * KX * KX;


				/* we are stepping through in row-major order */
				field_buffer[field_index(l, m, n, KX)] = field_val;
				
				// enforce real-ness of field
				//eprintf("%ld %ld %ld   " , l, m, n);
				field_buffer[field_index((KX - l) % KX, (KX - m) % KX, (KX - n) % KX, KX)]
					= conj(field_val);
				if ((l == 0 || (l == KX/2 && ! (KX % 2))) &&
						(m == 0 || (m == KX/2 && ! (KX % 2))) &&
						(n == 0 || (n == KX/2 && ! (KX % 2))) ){
					field_buffer[field_index(l, m, n, KX)] = creal(field_val);
				}
				#endif

			}
		}
	}

	/* Free the RNG to avoid a memory leak! */
	gsl_rng_free(rng);



}


