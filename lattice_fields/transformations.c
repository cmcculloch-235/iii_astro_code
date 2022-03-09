#include <math.h>
#include <complex.h>

#include "util.h"
#include "transformations.h"

void identity(void *in, void *out, size_t index, void *general_args)
{
	*(complex double *) out = *(complex double *) in;
}
