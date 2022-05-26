#ifndef MAT_H
#define MAT_H

double* affine_to_rot(double* A, int dim, bool hard_inv);
void sqmat_mult(double* A, double* B, double* AB, int dim);
void sqmat_t(double* A, double* At, int dim);
void eye(double* I, int dim);
void hard_inv_2d(double* A, double* Ai);

#endif