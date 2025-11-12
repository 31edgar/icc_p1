#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "linalg.h"

int main(void) {
    // Declaración de las variables
    int n;
    printf("Dame los elementos del sistema:");
    scanf("%d", &n);

    double tol, **A, **A_original, *b, *y, *r, *x, **L, **U;

    // Reservamos espacio en la memoria
    b = (double*)malloc(n*sizeof(double));
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

    L = (double**)malloc(n*sizeof(double*));
    for (int i=0; i<n; i++) 
        L[i] = (double*)malloc(n*sizeof(double));

    U = (double**)malloc(n*sizeof(double*));
    for (int i=0; i<n; i++) 
        U[i] = (double*)malloc(n*sizeof(double));
    
    // Nos aseguramos de que la memoria se ha asignado correctamente
    if (!x || !r || !b || !A || !U || !L) {
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

    // PEQUEÑOS AJUSTES EN LAS VARIABLES ANTES DE COMENZAR
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

    // Resolvemos el sistema utilizando la factorización LU
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

        }

        free(A_original);
        free(A);
        free(L);
        free(U);
        free(r);
        free(b);
        free(x);
        free(y);

        return 1;
    }

    // Si la podemos factorizar
    // Determinamos qué són L y U

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

    // Una vez tenemos L y U, solo queda resolver: Ly = b   y después   Ux = y
    
    // Ly = b
    // Este sistema lo podremos resolver siempre, porque los elementos de la diagonal de L son siempre 1
    resoltriinf(n, L, b, y);

    // Ahora hay que asegurarnos de que el otro sistema se puede resolver
    int errorResolucio = resoltrisup(n, U, y, x, tol); // 1 = error | 0 = correcto

    if (errorResolucio) {
        printf("ERROR: El sistema Ux = y no se pudo resolver.\n");

        // Liberamos la memoria usada
        for (int i=0; i<n; i++) {
            free(A[i]);
            free(L[i]);
            free(U[i]);
            free(A_original[i]);

        }

        free(A_original);
        free(A);
        free(L);
        free(U);
        free(r);
        free(b);
        free(x);
        free(y);

        return 1;
    }

    // Si hemos podido resolver los sistemas...
    // Calculamos el residuo que se guardará en r
    residuo(n, A_original, b, x, r);

    // Calculamos la norma2 del residuo
    double norma2 = sqrt(prod_esc(n, r, r));

    // Imprimimos la solución y la norma
    printf("Solución X:\n");
    for(int i=0; i<n; i++)
        printf("%lf ", x[i]);

    printf("\nNorma del residuo: %e\n", norma2);

    // Liberamos la memoria usada
        for (int i=0; i<n; i++) {
            free(A[i]);
            free(L[i]);
            free(U[i]);
            free(A_original[i]);

        }

        free(A);
        free(A_original);
        free(L);
        free(U);
        free(r);
        free(b);
        free(x);
        free(y);

    return 0;
}