#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdbool.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include "csparse.h"
#include "parser.h"
#include "list.h"
#include "hashtbl.h"
#include "table.h"
#include "util.h"

char *toString(int integer) {
    char *buffer = (char *) malloc(sizeof (integer));
    sprintf(buffer, "%d", integer);
    return buffer;
}

float calcItol(char *line) {
    char *token = strtok(line, " \n\r=");
    float itol = atof("1e-3");
    while (token != NULL) {
        if (strstr(token, "ITOL")) {
            itol = atof(strtok(NULL, " \n\r="));
            break;
        }
        token = strtok(NULL, " \n\r=");
    }
    return itol;
}

char *identifyElement(char *name) {
    if (name[0] == 'I' || name[0] == 'R' || name[0] == 'C' || name[0] == 'V' || name[0] == 'L'
            || name[0] == 'i' || name[0] == 'r' || name[0] == 'c' || name[0] == 'v' || name[0] == 'l')
        return "circuit";
    else if (strstr(name, "M") != NULL || strstr(name, "m") != NULL) return "mos";
    else if (strstr(name, "Q") != NULL || strstr(name, "q") != NULL) return "bjt";
    else if (strstr(name, ".DC") != NULL || strstr(name, ".dc") != NULL) return "dc";
    else if (strstr(name, ".TRAN") != NULL || strstr(name, ".tran") != NULL) return "tran";
    else if (strstr(name, ".PLOT") != NULL || strstr(name, ".plot") != NULL) return "plot";
}

int identifyCircuit(char *name) {
    if (strstr(name, "I") != NULL || strstr(name, "R") != NULL || strstr(name, "i") != NULL || strstr(name, "r") != NULL) return 1;
    else if (strstr(name, "V") != NULL || strstr(name, "v") != NULL) return 2;
}

void printVector(char *name, double *vector, size_t size) {
    printf("%s: \n", name);
    size_t i;
    for (i = 0; i < size; i++) {
        printf("%f\t\t", vector[i]);
    }
    printf("\n");
}

void printPlot(plotList *l, double *vector) {
    while (l != NULL) {
        printf("Plot: %f\n", vector[atoi(l->value) - 1]);
        printf("\n");
        l = l->nxt;
    }
}

float gslNorm(gsl_vector *vector) {
    size_t i;
    float sum = 0;
    for (i = 0; i < vector->size; i++) {
        sum += gsl_vector_get(vector, i) * gsl_vector_get(vector, i);
    }
    sum = sqrt(sum);
    return sum;
}

float norm(double *vector, size_t size) {
    size_t i;
    float sum = 0;
    for (i = 0; i < size; i++) {
        sum += vector[i] * vector[i];
    }
    sum = sqrt(sum);
    return sum;
}

double dotProduct(double *a, double *b, size_t size) {
    double sum = 0;
    size_t i;
    for (i = 0; i < size; i++) {
        sum += a[i] * b[i];
    }
    return sum;
}

double *mulConstant(double *a, double c, size_t size) {
    double *b = (double *) malloc(size * sizeof (double));
    size_t i;
    for (i = 0; i < size; i++) {
        b[i] = a[i] * c;
    }
    return b;
}

double *addVectors(double *a, double *b, size_t size) {
    double *c = (double *) malloc(size * sizeof (double));
    size_t i;
    for (i = 0; i < size; i++) {
        c[i] = a[i] + b[i];
    }
    return c;
}

cs *composeDiagonal(cs *A) {
    cs *M = cs_spalloc(A->m, A->n, A->nz, 1, 1);
    M->nz = A->nz;
    int i;
    int j = 0;
    for (i = 0; i < nonZeroElem; i++) {
        if (A->p[i] == A->i[i]) {
            M->p[j] = A->p[i];
            M->i[j] = A->i[i];
            if (!A->x[i]) M->x[j] = 1;
            else M->x[j] = A->x[i];
            //printf("i[%d]:%d\t", j, M->i[j]);
            //printf("p[%d]:%d\t", j, M->p[j]);
            //printf("x[%d]:%f\n", j, M->x[j]);
            j++;
        }

    }
    M->nz = j;
    return M;
}

double *solvePreconditioner(cs *M, double *z, double *r) {
    cs *M_inverse = cs_spalloc(M->m, M->n, M->nz, 1, 1);
    M_inverse->nz = M->nz;
    int i;
    for (i = 0; i < M_inverse->nz; i++) {
        M_inverse->p[i] = M->p[i];
        M_inverse->i[i] = M->i[i];
        M_inverse->x[i] = 1 / M->x[i];
        //printf("x[%d]:%d\t", i, M_inverse->x[i]);
        //printf("i[%d]:%d\t", i, M_inverse->i[i]);
        //printf("p[%d]:%d\t", i, M_inverse->p[i]);
    }
    M_inverse = cs_compress(M_inverse);
    cs_dupl(M_inverse);
    cs_gaxpy(M_inverse, r, z);
    //M_inverse = cs_compress(M_inverse);
    //cs_dupl(M_inverse);
    //cs_gaxpy(M_inverse, z, r);
    return z;
}

void printGSLMatrix(char *name, gsl_matrix *matrix) {
    printf("%s: \n", name);
    size_t i, j;
    for (i = 0; i < matrix->size1; i++) {
        for (j = 0; j < matrix->size2; j++) {
            printf("%.4g\t\t", gsl_matrix_get(matrix, i, j));
        }
        printf("\n");
    }
}

void printGslVector(char *name, gsl_vector *vector) {
    printf("%s: \n", name);
    size_t i;
    for (i = 0; i < vector->size; i++) {
        printf("%.4g\t\t", gsl_vector_get(vector, i));
    }
    printf("\n");
}

void printGslPlot(plotList *l, gsl_vector *vector) {
    while (l != NULL) {
        printf("Plot: %g\n", gsl_vector_get(vector, atoi(l->value) - 1));
        printf("\n");
        l = l->nxt;
    }
}

gsl_matrix *composeGslDiagonal(gsl_matrix *A) {
    gsl_matrix *M = gsl_matrix_calloc(A->size1, A->size2);
    int i;
    for (i = 0; i < M->size1; i++) {
        if (gsl_matrix_get(A, i, i) == 0) gsl_matrix_set(M, i, i, 1);
        else gsl_matrix_set(M, i, i, gsl_matrix_get(A, i, i));
    }
    return M;
}

gsl_vector *solveGslPreconditioner(gsl_matrix *M, gsl_vector *z, gsl_vector * r) {
    gsl_matrix *M_inverse = gsl_matrix_calloc(M->size1, M->size2);
    int i;
    for (i = 0; i < M_inverse->size1; i++) {
        if (gsl_matrix_get(M, i, i) != 0) gsl_matrix_set(M_inverse, i, i, 1 / gsl_matrix_get(M, i, i));
    }
    gsl_blas_dgemv(CblasNoTrans, 1.0, M_inverse, r, 0.0, z);
    //gsl_linalg_cholesky_decomp(M);
    //gsl_linalg_cholesky_solve(M, r, z);
    return z;
}

gsl_vector *computeGslAp(gsl_matrix *A, gsl_vector *p) {
    gsl_vector *q = gsl_vector_calloc(data + a2);
    gsl_blas_dgemv(CblasNoTrans, 1.0, A, p, 0.0, q);
    return q;
}

double *gslVectorToDouble(gsl_vector *gVector) {
    double *dVector = (double *) malloc((gVector->size) * sizeof (double));
    size_t i;
    for (i = 0; i < gVector->size; i++) {
        dVector[i] = gsl_vector_get(gVector, i);
    }
    return dVector;
}

gsl_vector *doubleToGslVector(double *dVector, size_t size) {
    gsl_vector *gVector = gsl_vector_calloc(size);
    size_t i;
    for (i = 0; i < size; i++) {
        gsl_vector_set(gVector, i, dVector[i]);
    }
    return gVector;
}

float sin_solve(float t, sinList * node) {
    //printf("sin_solve\n");
    //printf("t=%f\n", t);
    //printf("td=%f\n", atof(node->td));
    if (t >= 0 && t <= atof(node->td)) {
        //printf("sin_solve_if\n");
        return (atof(node->i1) + atof(node->ia) * sin(2 * (atof(node->ph) / 360) * M_PI));
    } else if (t > atof(node->td) && t <= atof(tran->end)) {
        return (atof(node->i1) + atof(node->ia) * sin(2 * atof(node->fr) * (t - atof(node->td)) * M_PI + 2 * (atof(node->ph) / 360) * M_PI) * exp(-(t - atof(node->td)) * atof(node->df)));
    }
}

float exp_solve(float t, expList * node) {
    if (t >= 0 && t <= atof(node->td1)) {
        return (atof(node->i1));
    } else if (t > atof(node->td1) && t <= atof(node->td2)) {
        return (atof(node->i1) + (atof(node->i2) - atof(node->i1))*(1 - exp(-(t - atof(node->td1)) / atof(node->tc1))));
    } else if (t > atof(node->td2) && t < atof(tran->end)) {
        return (atof(node->i1) + (atof(node->i2) - atof(node->i1))*(exp(-(t - atof(node->td2)) / atof(node->tc2)) - exp(-(t - atof(node->td1)) / atof(node->tc1))));
    }
}

float pulse_solve(int k, float t, pulseList * node) {
    if (t >= 0 && t <= atof(node->td)) {
        return (atof(node->i1));
    } else if (t >= atof(node->td) * k * atof(node->per) && t <= atof(node->td) + atof(node->tr) * k * atof(node->per)) {
        return (((atof(node->i1) - atof(node->i2)) / ((atof(node->td) * k * atof(node->per))-(atof(node->td) + atof(node->tr) * k * atof(node->per)))) * t);
    } else if (t >= atof(node->td) + atof(node->tr) * k * atof(node->per) && t <= atof(node->td) + atof(node->tr) + atof(node->pw) * k * atof(node->per)) {
        return (atof(node->i2));
    } else if (t >= atof(node->td) + atof(node->tr) + atof(node->pw) * k * atof(node->per) && t <= atof(node->td) + atof(node->tr) + atof(node->pw) + atof(node->tf) * k * atof(node->per)) {
        return ((atof(node->i2) - atof(node->i1)) / ((atof(node->td) + atof(node->tr) + atof(node->pw) + atof(node->tf) * k * atof(node->per))-(atof(node->td) + atof(node->tr) + atof(node->pw) * k * atof(node->per))) * t);
    } else if (t >= atof(node->td) + atof(node->tr) + atof(node->pw) + atof(node->tf) * k * atof(node->per) && t <= atof(node->td) + atof(node->per) * k * atof(node->per)) {
        return (atof(node->i1));
    }
}

float pwl_solve(float t, pwlList * node) {
    //node=pwlRoot;
    float t_start = 0;
    float t_end = atof(node->t1);
    //pwlList *cnt = node->pwl;
    while (node != NULL) {
        if (t >= t_start && t <= t_end) {
            return (((atof(node->i1)- atof(node->nxt->i1)) / (t_start-t_end)) * t);
        }
        node = node->nxt;
        t_start = t_end;
        t_end = atof(node->t1);
    }
}


void LU_solve(gsl_matrix *matrixA, gsl_vector *vectorX, gsl_vector *vectorB) {
    int s;
    gsl_permutation *p = gsl_permutation_alloc(data + a2);
    gsl_matrix *PA = gsl_matrix_calloc(data + a2, data + a2);
    gsl_matrix_memcpy(PA, matrixA);
    gsl_linalg_LU_decomp(PA, p, &s);
    gsl_linalg_LU_solve(PA, p, vectorB, vectorX);

}

