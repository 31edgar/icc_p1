#include <stdio.h>
#include <math.h>

#include "linalg.h"

int main(void){

    // Declaración de variables
    int n;
    scanf("%d", n);
    
    double *b, **A;

    b = (double*)malloc(n*sizeof(double));

    A = (double**)malloc(n*sizeof(double*));

    for (int i=0; i<n; i++) {
        A[i] = (double*)malloc(n*sizeof(double));
    }

    return 0;

    // Lectura de variables

    for(int i=0; i<n; i++){
        scanf("%lf", &b[i]);
    }

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            scanf("%lf", &A[i][j]);
        }
    }
}