/**
 *
 * Author: Marc Khoury <khoury@eecs.berkeley.edu>
 *
 * Copyright (C) 2014 Marc Khoury
 *
 * This version is released under the Eclipse Public License
 * with the Graphviz distribution.
 *
 * A version is also available on GitHub: https://github.com/marckhoury/mars
 * If you make improvements or bug fixes to this code it would be much
 * appreciated if you could also contribute those changes back to the
 * GitHub repository.
 */

#ifndef LINALG_H
#define LINALG_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define EPSILON 1E-6

typedef struct matrix {
    double* m;
    int r;
    int c;
}* mat;

mat mat_new(int r, int c);
mat mat_rand(int r, int c);
void mat_free(mat m);
void mat_print(mat m, FILE* out);
void mat_set(mat m, double val);
void mat_sub(mat a, mat b);
void mat_scalar_mult(mat m, double s);
mat mat_trans(mat m);
mat mat_mult(mat a, mat b);
double* mat_vec_mult(mat a, double* b);
double* mat_col(mat a, int j);
double mat_accu(mat m);
void vec_scalar_mult(double* v, int n, double s);
int mindex(int i, int j, mat m);
void vec_print(double* v, int n, FILE* out);
#endif
