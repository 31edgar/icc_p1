#include <stdio.h>
#include <math.h>

#include "linalg.h"

int main(void){

    int n;

    double *b, **A;

    b = (double*)malloc(n*sizeof(double));

    A = (double**)malloc(n*sizeof(double*));

    for (int i=0; i<n; i++) {
        A[i] = (double*)malloc(n*sizeof(double));
    }

    return 0;
}