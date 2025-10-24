#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "linalg.h"

int main(void){
    // Declaración de las variables
    int n;
    printf("Dame los elementos del sistema: ");
    scanf("%d", &n);

    double tol, *r, *b, *x, **A;

    // Reservamos espacio en la memoria para cada una de las variables
    x = (double*)malloc(n*sizeof(double));
    r = (double*)malloc(n*sizeof(double));
    b = (double*)malloc(n*sizeof(double));
    A = (double**)malloc(n*sizeof(double*));
    for (int i=0; i<n; i++)
        A[i] = (double*)malloc(n*sizeof(double));

    for (int i=0; i<n; i++)
        x[i] = 0;
    // Nos aseguramos de que la memoria se ha asignado correctamente
    if (!x || !r || !b || !A) {
        printf("Error al asignar memoria\n");
        return 1;
    }

    // Pedimos los datos
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

    // isTrisupResolta será 1 si hubo un fallo
    if (isTrisupResolta) {
        printf("La matriz no se pudo resolver como triangular superior.\n");
        // Liberamos la memoria antes de salir
        for(int i=0; i<n; i++) free(A[i]);
        free(A);
        free(b);
        free(x);
        free(r);
        return 1;
    }

    // Calculamos el residuo y este se guarda en r
    residuo(n, A, b, x, r);

    // Calculamos la norma2 del residuo
    double norma2 = sqrt(prod_esc(n, r, r));

    // Imprimimos la solución y la norma
    printf("Solución X:\n");
    for(int i=0; i<n; i++)
        printf("%lf ", x[i]);

    printf("\nNorma del residuo: %e\n", norma2);

    // Liberamos la memoria
    for(int i=0; i<n; i++) free(A[i]);
    free(A);
    free(b);
    free(x);
    free(r);

    return 0;
}
