#include <stdio.h>
#include <math.h>

#include "linalg.h"

int main(void){

    // Declaración de variables
    int n;
    printf("Dame los elementos del sistema");
    scanf("%d", n);

    double tol, *b, *x, **A;

    b = (double*)malloc(n*sizeof(double));

    A = (double**)malloc(n*sizeof(double*));
    for (int i=0; i<n; i++) {
        A[i] = (double*)malloc(n*sizeof(double));
    }

    return 0;

    // Lectura de variables
    printf("Dame las componentes del vector b\n");
    for(int i=0; i<n; i++){
        scanf("%lf", &b[i]);
    }

    printf("Dame las componentes de la matriz A\n");
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            scanf("%lf", &A[i][j]);
        }
    }

    tol = 1e-12;
}