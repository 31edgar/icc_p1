#include <stdio.h>
#include <math.h>
#include "linalg.h"

int resoltrisup(int n, double **a, double *b, double *x, double tol) {
    for (int i=n-1; i >= 0; i-- ) {
        double sum = 0.0;

        // Cálculo de la suma de todos los elementos a la derecha de la diagonal
        for (int k = i + 1; k < n; k++) {
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
    // Recorro las filas
    for (int i = 0; i < n; i++) {
        y[i] = 0.0; // Nos aseguramos que el vector donde guardamos el resultado és 0
        for (int j = 0; j < n; j++) {
            y[i] += a[i][j] * x[j]; // Calculamos cada elemento del vector
        }
    }
}

double prod_esc(int n, double *u, double *v) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += u[i] * v[i];
    }
    return sum;
}

void residuo(int n, double **a, double *b, double *x, double *r) {
    // Guardem espai a memoria de forma dinàmica i ens assegurem de que sigui possible fer-ho
    double *Ax = (double *)malloc(n * sizeof(double));
    if (Ax == NULL) return;

    // Multipliquem la nostra matriu pel vector
    prodMatVec(n, a, x, Ax);  // Ax = A * x

    // Guardem al vector residu les diferencies entre els productes y els termes independents (RESIDU)
    for (int i = 0; i < n; i++) {
        r[i] = Ax[i] - b[i];  // r = A*x - b
    }

    // Cal recordar alliberar la memòria que hem ocupat a la primera línia de la funció
    free(Ax);
}



