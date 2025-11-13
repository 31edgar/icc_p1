#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "linalg.h"

double resolverSistema(int n, double **A, double *b, double *r, double *x, double tol) {
    // Copiamos A
    double **A_original = (double**)malloc(n*sizeof(double*));
    for (int i=0; i<n; i++)
        A_original[i] = (double*)malloc(n*sizeof(double));

    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++)
            A_original[i][j] = A[i][j];
    }

    // Resolvemos el sistema
    int errorEliminacio = gausspiv(n, A, b, tol);

    // isTrisupResolta será 1 si hubo un fallo
    if (errorEliminacio) {
        printf("ERROR: no se pudo aplicar el método de eliminación por Gauss\n");
        return -1;
    }

    int errorResolució = resoltrisup(n, A, b, x, tol);

    if (errorResolució) {
        printf("ERROR: no se pudo resolver el sistema. \n");
        // Liberamos la memoria antes de salir
        return -1;
    }

    // Calculamos el residuo y este se guarda en r
    residuo(n, A_original, b, x, r);

    // Calculamos la norma2 del residuo
    double norma2 = sqrt(prod_esc(n, r, r));

    return norma2;
}

int main(void) {
    srand(time(NULL));
    double tol = 1e-12;
    double norma2max, norma2;

    // Preparamos el archivo donde vamos a escribir los datos
    FILE *f = fopen("salida.txt", "w");
    if (f == NULL) {
        perror("No se pudo abrir el archivo");
        return 1;
    }

    fprintf(f, "Dimensió    Norma 2 màxima\n");

    // Las 100 dimensiones diferentes
    for (int n=1; n<101; n++) {
        // Reservamos espacio para la matriz que vamos a utilizar
        double *r = (double*)malloc(n*sizeof(double));
        double *x = (double*)malloc(n*sizeof(double));
        double *b = (double*)malloc(n*sizeof(double));
        double **A = (double**)malloc(n*sizeof(double*));
        for (int k=0; k<n; k++)
            A[k] = (double*)malloc(n*sizeof(double));

        norma2max = 0.0;

        // Por cada dimension, resolvemos 10 sistemas
        for (int k=0; k<10; k++) {
            // Escribimos los valores de A
            for (int i=0; i<n; i++)
                for (int j=0; j<n; j++)
                    // El número debe ser entre [-1, 1]
                    A[i][j] = -1.0 + 2.0 * ((double)rand() / RAND_MAX);

            // Escribimos los valores de B
            for (int i=0; i<n; i++)
                b[i] = -1.0 + 2.0 * ((double)rand() / RAND_MAX);

            norma2 = resolverSistema(n, A, b, r, x, tol);
            if (norma2 > norma2max) norma2max = norma2;
        }

        fprintf(f, "%d %f\n", n, norma2max);


        for (int k=0; k<n; k++)
            free(A[k]);
        free(A);
        free(b);
        free(x);
        free(r);
    }

    fclose(f);

    return 0;
}