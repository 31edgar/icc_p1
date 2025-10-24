#ifndef LINALG_H
#define LINALG_H

#include <stdio.h>
#include <math.h>

// Declaración de funciones
int resoltrisup(int n, double **a, double *b, double *x, double tol);
void residuo(int n, double **a, double *b, double *x, double *r);
void prodMatVec(int n, double **a, double *x, double *y);
double prod_esc(int n, double *u, double *v);
void resoltriinf(int n, double **a, double *b, double *x);


#endif