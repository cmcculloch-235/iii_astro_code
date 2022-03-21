#include <math.h>
#include <complex.h>

#include "thread_pool.h"
#include "util.h"
#include "transformations.h"

void transform_identity(complex double *field, size_t N, size_t N_THREADS )
{
	return;
}

static void exp_map(void *in, void *out, size_t index, void *gen_args);
static void power_depth_map(void *in, void *out, size_t index, void *gen_args);
static void contrast(complex double *field, size_t N, size_t N_THREADS);

void transform_log_normal_contrast(complex double *field, size_t N, size_t N_THREADS)
{
	/* See notes for definition of this model */
	

	thread_map(&exp_map, (void *) NULL,
			(void *) field, sizeof(complex double),
			(void *) field, sizeof(complex double),
			N, N_THREADS);

	contrast(field, N, N_THREADS);
	
	thread_map(&power_depth_map, (void *) NULL,
			(void *) field, sizeof(complex double),
			(void *) field, sizeof(complex double),
			N, N_THREADS);

	contrast(field, N, N_THREADS);

}

void transform_power_law_contrast(complex double *field, size_t N, size_t N_THREADS)
{

	thread_map(&power_depth_map, (void *) NULL,
			(void *) field, sizeof(complex double),
			(void *) field, sizeof(complex double),
			N, N_THREADS);

	contrast(field, N, N_THREADS);
}


static void exp_map(void *in, void *out, size_t index, void *gen_args)
{

	*(complex double *) out = cexp(*(complex double *) in);

}

static const double A = 1.0;
static const double ALPHA = 1.7;
static void power_depth_map(void *in, void *out, size_t index, void *gen_args)
{
	double in_val = creal(*(complex double *) in);
	if (in_val < 0.0) {
		in_val = 0.0;
	}
	*(complex double *) out = A * cpow(1.0 + in_val, ALPHA);

}

static void divide_by_mean_and_sub1(void *in, void *out, size_t index, void *general_args);

static void contrast(complex double *field, size_t N, size_t N_THREADS)
{
	/* Transform the field in place to its contrast */

	/* find the mean value without parallelisation as that's the easiest
	 * way to be thread-safe */
	complex double field_mean = 0.0;
	for (size_t i = 0; i < N; ++i) {
		field_mean += field[i];
	}
	field_mean /= N;

	thread_map(&divide_by_mean_and_sub1, (void *) &field_mean,
			(void *) field, sizeof(complex double),
			(void *) field, sizeof(complex double),
			N, N_THREADS);
}

static void divide_by_mean_and_sub1(void *in, void *out, size_t index, void *mean_ptr)
{
	complex double mean = *(complex double *) mean_ptr;
	*(complex double *) out = (*(complex double *) in) / mean - 1.0;
}


