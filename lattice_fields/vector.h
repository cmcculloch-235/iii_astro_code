#ifndef VEC_H_INC
#define VEC_H_INC

/* I wrote this years ago for another thing, so it has some extra features */
void vec_print(double a[3]);

double vec_dot(double a[3], double b[3]);
void vec_sca_mul(double a[3], double k, double *out);
void vec_add(double a[3], double b[3], double *out);
void vec_sub(double a[3], double b[3], double *out);
void vec_cross(double a[3], double b[3], double *out);
void vec_rotate(double a[3], double axis_nu[3], double theta,
		double *out);

#endif
