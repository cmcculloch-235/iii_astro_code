#include <math.h>
#include <stdio.h>
#include "vector.h"

// The purpose of this is to work with vectors of the double type.

void vec_print(double a[3]) 
{
	double mag = pow(vec_dot(a, a), 0.5);
	printf("(%f, %f, %f), %f\n", a[0], a[1], a[2], mag);
}

double vec_dot(double a[3], double b[3])
{
	double r = 0;
	for (int i = 0; i < 3; ++i) {
		r += a[i] * b[i];
	}
	return r;
}

// Multiplies the components of the vector a by the scalar k, storing the
// result in the vector out
void vec_sca_mul(double a[3], double k, double *out)
{
	for (int i = 0; i < 3; ++i) {
		out[i] = a[i] * k;
	}
}

// Adds the vector a to the vector b, storing the result in the vector out
void vec_add(double a[3], double b[3], double *out)
{
	for (int i = 0; i < 3; ++i) {
		out[i] = a[i] + b[i];
	}
}

/* Vector subtraction */
void vec_sub(double a[3], double b[3], double *out)
{
	vec_sca_mul(b, -1, out);
	vec_add(a, out, out);
}

// Computes the vector product of the vectors a and b, storing the result in
// the vector out. Has a temporary vector so that if *out points to a or b,
// unexpected behaviour does not result.
void vec_cross(double a[3], double b[3], double *out)
{
	double tmp[3];
	for (int i = 0; i < 3; ++i) {
		tmp[i] = a[(i + 1) % 3] * b[(i + 2) % 3]
			- a[(i + 2) % 3] * b[(i + 1) % 3];
	}
	for (int i = 0; i < 3; ++i) {
		out[i] = tmp[i];
	}
}

void vec_rotate(double a[3], double axis_nu[3], double theta, double *out)
{
	// As discussed in analysis.

	double axis[3];
	// Find a unit vector in the direction of axis_nu
	vec_sca_mul(axis_nu, 1/pow(vec_dot(axis_nu, axis_nu), 0.5), axis);
	// Find the projection of a onto axis
	double beta = vec_dot(a, axis);
	
	// Find the component of vector a perpendicular to the axis
	double a_naxis[3];
	double tmp[3][3];
	vec_sca_mul(axis, -1 * beta, tmp[0]);
	vec_add(a, tmp[0], a_naxis);
	double alpha = pow(vec_dot(a_naxis, a_naxis), 0.5);
	
	// Find a unit vector in the direction of a perpendicular to the axis
	double p_hat[3];
	vec_sca_mul(a_naxis, 1/alpha, p_hat);
	// Find the cross product of the axis and the unit vector perpendicular to
	// it
	double axis_x_p_hat[3];
	vec_cross(axis, p_hat, axis_x_p_hat);
	
	// a' = beta*axis + alpha * (p_hat cos theta + (axis x p_hat) sin theta)
	vec_sca_mul(axis, beta, tmp[0]);
	vec_sca_mul(p_hat, alpha * cos(theta), tmp[1]);
	vec_sca_mul(axis_x_p_hat, alpha * sin(theta), tmp[2]);

	for (int i = 0; i < 3; ++i) {
		out[i] = tmp[0][i] + tmp[1][i] + tmp[2][i];
	}
}
