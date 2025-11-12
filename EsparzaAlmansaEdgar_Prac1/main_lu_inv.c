#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "linalg.h"

int main(void) {
    // Declaración de las variables
    int n;
    printf("Dame los elementos del sistema:");
    scanf("%d", &n);

    double tol, **A, **A_original, **A_inv, *y, *r, *x, **L, **U;

    // Reservamos espacio en la memoria
    y = (double*)malloc(n*sizeof(double));
    r = (double*)malloc(n*sizeof(double));

    x = (double*)malloc(n*sizeof(double));
    for (int i=0; i<n; i++)
        x[i] = 0;

    A = (double**)malloc(n*sizeof(double*));
    for (int i=0; i<n; i++) 
        A[i] = (double*)malloc(n*sizeof(double));

    A_original = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++)
        A_original[i] = (double*)malloc(n*sizeof(double));

    A_inv = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++)
        A_inv[i] = (double*)malloc(n*sizeof(double));

    L = (double**)malloc(n*sizeof(double*));
    for (int i=0; i<n; i++) 
        L[i] = (double*)malloc(n*sizeof(double));

    U = (double**)malloc(n*sizeof(double*));
    for (int i=0; i<n; i++) 
        U[i] = (double*)malloc(n*sizeof(double));
    
    // Nos aseguramos de que la memoria se ha asignado correctamente
    if (!x || !r || !A || !U || !L || !A_original || !A_inv) {
        printf("Error al asignar memoria\n");
        return 1;
    }

    // Pedimos los datos
    printf("Dame las componentes de la matriz A\n");
    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++)
            scanf("%lf", &A[i][j]);
    
    printf("Introduce la tolerancia: ");
    scanf("%lf", &tol);



    // PEQUEÑOS AJUSTES EN LAS VARIABLES ANTES DE COMENZAR ----------------------------------
    // Guardamos la A original ya que cuando llamemos a la función lu() cambiaremos sus valores
    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++)
            A_original[i][j] = A[i][j];

    // Inicializamos L y U en 0
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            U[i][j] = 0;
            L[i][j] = 0;
        }
    }      

    // INVERSIÓN DE LA MATRIZ --------------------------------------------------------------
    // Factorizamos usando LU
    int errorFactoritzacio = lu(n, A, tol); // 1 = error | 0 = correcto

    // CASO EN EL QUE LA MATRIZ NO SE PUEDE FACTORIZAR
    if (errorFactoritzacio) {
        printf("ERROR: La matriz no se ha podido factorizar.\n");

        // Liberamos la memoria usada
        for (int i=0; i<n; i++) {
            free(A[i]);
            free(L[i]);
            free(U[i]);
            free(A_original[i]);
            free(A_inv[i]);
        }

        free(A);
        free(A_original);
        free(L);
        free(U);
        free(r);
        free(A_inv);
        free(x);
        free(y);

        return 1;
    }

    // Extraemos L y U de A (A ahora es una "fusión" de las dos)
    // U son los elementos de la diagonal y lo que hay arriba
    for (int i=0; i<n; i++) {
        for (int j=i; j<n; j++) {
            U[i][j] = A[i][j];
        }
    }

    // L son los elementos que hay abajo de la diagonal (L tiene una diagonal de 1s)
    // Si nos fijamos, recorremos la matriz igual que para U pero intercambiando I por J
    for (int j=0; j<n; j++) {
        for (int i=j; i<n; i++) {
            if (i == j) {
                L[i][j] = 1; // Si es una elemento de la diagonal, no lo copiaremos y pondremos 1
            } else {
                L[i][j] = A[i][j]; // Copiamos el elemento
            }
        }
    }

    // Una vez tenemos L y U, solo queda resolver: Ly = e(i)   y   Ux(i) = y   unas n veces
    // En cada iteración, conseguiremos x que es un vector, este vector es nada más ni nada menos que
    // una columna de la matriz inversa, en este caso la columna i

    // Para cada columna i de la inversa:
    for (int i = 0; i < n; i++) {
        double *e = (double*)calloc(n, sizeof(double)); // vector e_i (todo 0s)
        
        e[i] = 1; // el vector base canónico correspondiente

        // Resolver Ly = e(i)
        resoltriinf(n, L, e, y);

        // Resolver Ux(i) = y
        int errorResolucio = resoltrisup(n, U, y, x, tol);

        if (errorResolucio) {
            printf("ERROR: No se ha podido solucionar el sistema en la iteración %d\n", i);

            // Liberamos la memoria usada
            for (int i=0; i<n; i++) {
                free(A[i]);
                free(L[i]);
                free(U[i]);
                free(A_original[i]);
                free(A_inv[i]);
            }

            free(A);
            free(A_original);
            free(L);
            free(U);
            free(r);
            free(A_inv);
            free(x);
            free(y);

            return 1;
        }

        // Guardar x en la columna i de A_inv
        for (int j = 0; j < n; j++) {
            A_inv[j][i] = x[j];
        }

        free(e);
    }

    // ------------------------------------------------------------------------------------
    // CÁLCULO DE AA⁻1 - Id

    A = prodMatMat(n, n, n, A_original, A_inv);

    for (int i=0; i<n; i++)
        A[i][i] -= 1;

    // CÁLCULO DE LA NORMA INFINITO -------------------------------------------------------
    // Queremos calcular: la suma de los valores absolutos de cada fila

    double maxSum = tol; // Empezamos con un número muy pequeño, por ejemplo, la tolerancia

    for (int i=0; i<n; i++) {
        double sum = 0;

        for (int j=0; j<n; j++)
            sum += fabs(A[i][j]);

        if (sum > maxSum) maxSum = sum;
    }


    // ------------------------------------------------------------------------------------
    // Imprimimos A⁻1
    printf("A inversa:\n");

    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++)
            printf("%lf ", A_inv[i][j]);
        printf("\n");
    }

    printf("Norma infinita: %lf", maxSum);

    // Liberamos la memoria usada
    for (int i=0; i<n; i++) {
        free(A[i]);
        free(L[i]);
        free(U[i]);
        free(A_original[i]);
        free(A_inv[i]);
    }

    free(A);
    free(A_original);
    free(L);
    free(U);
    free(r);
    free(A_inv);
    free(x);
    free(y);

    return 0;

    
}