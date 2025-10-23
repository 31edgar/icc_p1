#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "linalg.h"

int main(void){
    int n;
    printf("Dame los elementos del sistema: ");
    scanf("%d", &n);

    double tol, *r, *b, *x, **A;

    x = (double*)malloc(n*sizeof(double));
    r = (double*)malloc(n*sizeof(double));
    b = (double*)malloc(n*sizeof(double));
    A = (double**)malloc(n*sizeof(double*));
    for (int i=0; i<n; i++)
        A[i] = (double*)malloc(n*sizeof(double));

    if (!x || !r || !b || !A) {
        printf("Error al asignar memoria\n");
        return 1;
    }

    printf("Dame las componentes del vector b\n");
    for(int i=0; i<n; i++)
        scanf("%lf", &b[i]);

    printf("Dame las componentes de la matriz A\n");
    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++)
            scanf("%lf", &A[i][j]);

    printf("Introduce la tolerancia: ");
    scanf("%lf", &tol);

    int isTrisupResolta = resoltrisup(n, A, b, x, tol);
    if (isTrisupResolta) {  // si es 1, hubo un fallo
        printf("La matriz no se pudo resolver como triangular superior.\n");
        // liberar memoria antes de salir
        for(int i=0; i<n; i++) free(A[i]);
        free(A);
        free(b);
        free(x);
        free(r);
        return 1;
    }

    residuo(n, A, b, x, r);
    double norma2 = sqrt(prod_esc(n, r, r));

    printf("Solución X:\n");
    for(int i=0; i<n; i++)
        printf("%lf ", x[i]);
    printf("\nNorma del residuo: %e\n", norma2);

    for(int i=0; i<n; i++) free(A[i]);
    free(A);
    free(b);
    free(x);
    free(r);

    return 0;
}
