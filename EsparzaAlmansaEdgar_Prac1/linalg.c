#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "linalg.h"

// 1. Resolución de sistemas triangulares
int resoltrisup(int n, double **a, double *b, double *x, double tol) {
    int i, k;

    for (i=n-1; i >= 0; i-- ) {
        double sum = 0.0;

        // Cálculo de la suma de todos los elementos a la derecha de la diagonal
        for (k = i + 1; k < n; k++) {
            sum += a[i][k] * x[k];
        }

        // Comprobamos que el número por el que dividiremos sea suficientemente grande
        if (fabs(a[i][i]) < tol)
            return 1; // No se puede dividir por algo casi 0

        // Calculamos y guardamos en el vector x el valor de la incógnita
        x[i] = (b[i] - sum) / a[i][i];
    }

    return 0;
}

void prodMatVec(int n, double **a, double *x, double *y) {
    int i, j;

    // Recorro las filas
    for (i = 0; i < n; i++) {
        y[i] = 0.0; // Nos aseguramos que el vector donde guardamos el resultado és 0
        for (j = 0; j < n; j++) {
            y[i] += a[i][j] * x[j]; // Calculamos cada elemento del vector
        }
    }
}

double **prodMatMat(int n, int p, int m, double **A, double **B) {
    // Reservem espai per la matriu resultat C (n x m)
    double **C = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++)
        C[i] = (double*)malloc(m * sizeof(double));

    // Calculem C = A * B
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            C[i][j] = 0.0;
            for (int k = 0; k < p; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    // Missatge informatiu sobre si BA es pot calcular
    if (m == n)
        printf("També és possible calcular BA (donarà una matriu de %dx%d)\n", p, p);
    else
        printf("No és possible calcular BA, perquè el nombre de columnes de B (%d) no coincideix amb les files d'A (%d)\n", m, n);

    return C;
}

double prod_esc(int n, double *u, double *v) {
    int i;
    double sum = 0.0;

    // Multiplicamos el componente i de un vector con el componente i del otro vector
    for (i = 0; i < n; i++) {
        sum += u[i] * v[i];
    }

    return sum;
}

void residuo(int n, double **a, double *b, double *x, double *r) {
    int i;

    // Guardem espai a memoria de forma dinàmica i ens assegurem de que sigui possible fer-ho
    double *Ax = (double *)malloc(n * sizeof(double));
    if (Ax == NULL) return;

    // Multipliquem la nostra matriu pel vector
    prodMatVec(n, a, x, Ax);  // Ax = A * x

    // Guardem al vector residu les diferencies entre els productes y els termes independents (RESIDU)
    for (i = 0; i < n; i++) {
        r[i] = Ax[i] - b[i];  // r = A*x - b
    }

    // Cal recordar alliberar la memòria que hem ocupat a la primera línia de la funció
    free(Ax);
}

void resoltriinf(int n, double **a, double *b, double *x) {
    // Asumimos que utilizaremos esta función para resolver una matriz triangular inferior con 1s en la diagonal
    int i, j;

    // Recorremos cada fila
    for (i=0; i<n; i++) {
        double sum = 0.0;

        // Sumamos todo lo de antes de la diagonal (multiplicado por su respectiva solucion Xj)
        for (j=0; j<i; j++) {
            sum += a[i][j] * x[j];
        }

        // Una vez tenemos la suma anterior, restamos ese valor al término independiente de la fila
        x[i] = b[i] - sum;
    }
}

// 2. Factorización LU
int lu(int n, double **a, double tol) {
    int i, j, k;
    double mult;

    // K representa la fila/columna donde se encuentra el pivote
    for (k = 0; k < n - 1; k++) {
        // Comprobamos que el pivote no sea demasiado pequeño
        if (fabs(a[k][k]) < tol) {
            return 1; // En ese caso no se puede continuar
        }

        // i representa cada fila de debajo del pivote
        for (i = k + 1; i < n; i++) {
            mult = a[i][k] / a[k][k]; // Calculamos el multiplicador
            a[i][k] = mult; // Lo almacenamos en la parte inferior (L)

            // J representa cada columna a la derecha del pivote
            for (j = k + 1; j < n; j++) {
                a[i][j] -= mult * a[k][j];
            }
        }

        // En resumen, para cada pivote[k][k], iteramos por todas sus filas i
        // A cada una, calculamos su multplicador (mult) y para cada elemento j de la fila i
        // Le restamos el producto (mult*elemento j de la fila k, fila anterior a i) (k=i-1)
    }

    return 0; // Éxito
}

// 3. Método de Gauss con pivoteo máximo por columnas
int gausspiv(int n, double **a, double *b, double tol) {
    for (int k = 0; k < n - 1; k++) {
        // K es la fila donde debe estar el pivote
        // i es la fila donde encontraremos el número con valor absoluto más alto (verdadero pivote)
        int pivoteFila = k;
        double maxVal = fabs(a[k][k]);
        for (int i = k + 1; i < n; i++) {
            if (fabs(a[i][k]) > maxVal) {
                maxVal = fabs(a[i][k]);
                pivoteFila = i;
            }
        }

        // Comprobar si el pivote es menor que la tolerancia
        if (maxVal < tol) {
            // Matriz singular o muy mal condicionada
            return 0; // Error, no se puede continuar
        }

        // Intercambiar filas si pivote máximo no está en la fila actual k
        if (pivoteFila != k) {
            // Intercambio de filas en la matriz A
            double *tempFila = a[k];
            a[k] = a[pivoteFila];
            a[pivoteFila] = tempFila;

            // Intercambio en el vector b
            double tempB = b[k];
            b[k] = b[pivoteFila];
            b[pivoteFila] = tempB;
        }

        // Eliminación hacia abajo
        for (int i = k + 1; i < n; i++) {
            double factor = a[i][k] / a[k][k];
            a[i][k] = 0.0; // explícito para claridad
            for (int j = k + 1; j < n; j++) {
                a[i][j] -= factor * a[k][j];
            }
            b[i] -= factor * b[k];
        }
    }

    // Si se llega aquí, eliminación gaussiana con pivoteo máximo fue exitosa
    return 1;
}