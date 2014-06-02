#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdbool.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include "csparse.h"
#include "parser.h"
#include "list.h"
#include "hashtbl.h"
#include "table.h"
#include "util.h"


gsl_matrix *A;
gsl_matrix *diag_A;
gsl_matrix *tilda_G;
gsl_matrix *tilda_C;
gsl_vector *x;
gsl_vector *b;
gsl_vector *tran_x;
gsl_vector *tran_b;
gsl_vector *tran_b_nxt;

cs *cs_A;
cs *cs_C;
double *cs_b;
double *cs_x;

HASHTBL *hashtbl;

int data = 0;
int nonZeroElem = 0;
int a1 = 0;
int a2 = 0;

HASHTBL *populateHashtbl(circuitList *node) {
    //printf("*** populateHashtbl ***\n");
    HASHTBL *hashtbl = hashtbl_create(65536, NULL);
    while (node != NULL) {
        if (identifyCircuit(node->name) == 1) a1++;
        else if (identifyCircuit(node->name) == 2) a2++;
        if (!strcmp(node->positive, "0")) {
            hashtbl_insert(hashtbl, node->positive, "0");
        } else if (hashtbl_insert(hashtbl, node->positive, toString(data + 1))) {
            data++;
        }
        if (!strcmp(node->negative, "0")) {
            hashtbl_insert(hashtbl, node->negative, "0");
        } else if (hashtbl_insert(hashtbl, node->negative, toString(data + 1))) {
            data++;
        }
        node = node->nxt;
    }
    //printf("\n");
    return hashtbl;
}

void allocateSystem() {
    A = gsl_matrix_calloc(data + a2, data + a2);
    x = gsl_vector_alloc(data + a2);
    b = gsl_vector_calloc(data + a2);
}

void allocateCSSystem() {
    cs_A = cs_spalloc(data + a2, data + a2, nonZeroElem, 1, 1);
    cs_A->nz = nonZeroElem;
    if (todoCSCG || todoCSBiCG) diag_A = gsl_matrix_calloc(data + a2, data + a2);
    cs_b = (double *) malloc((data + a2) * sizeof (double));
    cs_x = (double *) malloc((data + a2) * sizeof (double));
}

void allocateTransientSystem() {
    //printf("alloc trans\n");
    tilda_G = gsl_matrix_calloc(data + a2, data + a2);
    tilda_C = gsl_matrix_calloc(data + a2, data + a2);
    tran_x = gsl_vector_calloc(data + a2);
    tran_b = gsl_vector_calloc(data + a2);
    tran_b_nxt = gsl_vector_calloc(data + a2);
}

void populateSystem(circuitList *n, HASHTBL *hashtbl, gsl_matrix * A, gsl_vector * b) {
    int positive;
    int negative;
    float prevVal;
    int Vcnt = 0;
    while (n != NULL) {
        positive = atoi(hashtbl_get(hashtbl, n->positive));
        negative = atoi(hashtbl_get(hashtbl, n->negative));
        if (strstr(n->name, "R") || strstr(n->name, "r")) {
            if (!positive) {
                prevVal = gsl_matrix_get(A, negative - 1, negative - 1);
                gsl_matrix_set(A, negative - 1, negative - 1, 1 / atof(n->value) + prevVal);
            } else if (!negative) {
                prevVal = gsl_matrix_get(A, positive - 1, positive - 1);
                gsl_matrix_set(A, positive - 1, positive - 1, 1 / atof(n->value) + prevVal);
            } else {
                prevVal = gsl_matrix_get(A, positive - 1, positive - 1);
                gsl_matrix_set(A, positive - 1, positive - 1, 1 / atof(n->value) + prevVal);
                prevVal = gsl_matrix_get(A, positive - 1, negative - 1);
                gsl_matrix_set(A, positive - 1, negative - 1, -(1 / atof(n->value)) + prevVal);
                prevVal = gsl_matrix_get(A, negative - 1, positive - 1);
                gsl_matrix_set(A, negative - 1, positive - 1, -(1 / atof(n->value)) + prevVal);
                prevVal = gsl_matrix_get(A, negative - 1, negative - 1);
                gsl_matrix_set(A, negative - 1, negative - 1, 1 / atof(n->value) + prevVal);
            }
        } else if (strstr(n->name, "I") || strstr(n->name, "i")) {
            if (!positive) {
                prevVal = gsl_vector_get(b, negative - 1);
                gsl_vector_set(b, negative - 1, atof(n->value) + prevVal);
            } else if (!negative) {
                prevVal = gsl_vector_get(b, positive - 1);
                gsl_vector_set(b, positive - 1, -atof(n->value) + prevVal);
            } else {
                prevVal = gsl_vector_get(b, positive - 1);
                gsl_vector_set(b, positive - 1, -atof(n->value) + prevVal);
                prevVal = gsl_vector_get(b, negative - 1);
                gsl_vector_set(b, negative - 1, atof(n->value) + prevVal);
            }
        } else if (strstr(n->name, "V") || strstr(n->name, "v")) {
            if (!positive) {
                gsl_matrix_set(A, data + Vcnt, negative - 1, -1);
                gsl_matrix_set(A, negative - 1, data + Vcnt, -1);
            } else if (!negative) {
                gsl_matrix_set(A, data + Vcnt, positive - 1, 1);
                gsl_matrix_set(A, positive - 1, data + Vcnt, 1);
            } else {
                gsl_matrix_set(A, data + Vcnt, positive - 1, 1);
                gsl_matrix_set(A, data + Vcnt, negative - 1, -1);
                gsl_matrix_set(A, positive - 1, data + Vcnt, 1);
                gsl_matrix_set(A, negative - 1, data + Vcnt, -1);
            }
            gsl_vector_set(b, data + Vcnt, atof(n->value));
            Vcnt++;
        }
        n = n->nxt;
    }
    //printGSLMatrix("A", A);
    //printf("\n");
    //printGslVector("b", b);
}

void populateCSSystem(circuitList *n, HASHTBL *hashtbl, cs * cs_A) {
    int positive, negative;
    int Vcnt = 0;
    float prevVal;
    int i = 0;
    int j;
    //printf("\n");
    while (n != NULL) {
        positive = atoi(hashtbl_get(hashtbl, n->positive));
        negative = atoi(hashtbl_get(hashtbl, n->negative));
        if (strstr(n->name, "R") || strstr(n->name, "r")) {
            if (!positive) {
                cs_A->i[i] = negative - 1;
                cs_A->p[i] = negative - 1;
                cs_A->x[i] = 1 / atof(n->value);
                //printf("[%d, %d]: %.4f\n", cs_A->i[i], cs_A->p[i], cs_A->x[i]);
                i++;
                if (todoCSCG || todoCSBiCG) {
                    prevVal = gsl_matrix_get(diag_A, negative - 1, negative - 1);
                    if (1 / atof(n->value) + prevVal != 0) gsl_matrix_set(diag_A, negative - 1, negative - 1, 1 / atof(n->value) + prevVal);
                    else gsl_matrix_set(diag_A, negative - 1, negative - 1, 1);
                }
            } else if (!negative) {
                cs_A->i[i] = positive - 1;
                cs_A->p[i] = positive - 1;
                cs_A->x[i] = 1 / atof(n->value);
                //printf("[%d, %d]: %.4f\n", cs_A->i[i], cs_A->p[i], cs_A->x[i]);
                i++;
                if (todoCSCG || todoCSBiCG) {
                    prevVal = gsl_matrix_get(diag_A, positive - 1, positive - 1);
                    if (1 / atof(n->value) + prevVal != 0) gsl_matrix_set(diag_A, positive - 1, positive - 1, 1 / atof(n->value) + prevVal);
                    else gsl_matrix_set(diag_A, positive - 1, positive - 1, 1);
                }
            } else {
                cs_A->i[i] = positive - 1;
                cs_A->p[i] = positive - 1;
                cs_A->x[i] = 1 / atof(n->value);
                //printf("[%d, %d]: %.4f\n", cs_A->i[i], cs_A->p[i], cs_A->x[i]);
                i++;
                cs_A->i[i] = positive - 1;
                cs_A->p[i] = negative - 1;
                cs_A->x[i] = -1 / atof(n->value);
                //printf("[%d, %d]: %.4f\n", cs_A->i[i], cs_A->p[i], cs_A->x[i]);
                i++;
                cs_A->i[i] = negative - 1;
                cs_A->p[i] = positive - 1;
                cs_A->x[i] = -1 / atof(n->value);
                //printf("[%d, %d]: %.4f\n", cs_A->i[i], cs_A->p[i], cs_A->x[i]);
                i++;
                cs_A->i[i] = negative - 1;
                cs_A->p[i] = negative - 1;
                cs_A->x[i] = 1 / atof(n->value);
                //printf("[%d, %d]: %.4f\n", cs_A->i[i], cs_A->p[i], cs_A->x[i]);
                i++;
                if (todoCSCG || todoCSBiCG) {
                    prevVal = gsl_matrix_get(diag_A, positive - 1, positive - 1);
                    if (1 / atof(n->value) + prevVal != 0) gsl_matrix_set(diag_A, positive - 1, positive - 1, 1 / atof(n->value) + prevVal);
                    else gsl_matrix_set(diag_A, positive - 1, positive - 1, 1);
                    prevVal = gsl_matrix_get(diag_A, negative - 1, negative - 1);
                    if (1 / atof(n->value) + prevVal != 0) gsl_matrix_set(diag_A, negative - 1, negative - 1, 1 / atof(n->value) + prevVal);
                    else gsl_matrix_set(diag_A, negative - 1, negative - 1, 1);
                }
            }
        } else if (strstr(n->name, "I") || strstr(n->name, "i")) {
            if (!positive) {
                prevVal = cs_b[negative - 1];
                cs_b[negative - 1] = atof(n->value) + prevVal;
            } else if (!negative) {
                prevVal = cs_b[positive - 1];
                cs_b[positive - 1] = -atof(n->value) + prevVal;
            } else {
                prevVal = cs_b[positive - 1];
                cs_b[positive - 1] = -atof(n->value) + prevVal;
                prevVal = cs_b[negative - 1];
                cs_b[negative - 1] = atof(n->value) + prevVal;
            }
        } else if (strstr(n->name, "V") || strstr(n->name, "v")) {
            if (!positive) {
                cs_A->i[i] = data + Vcnt;
                cs_A->p[i] = negative - 1;
                cs_A->x[i] = -1;
                //printf("[%d, %d]: %.4f\n", cs_A->i[i], cs_A->p[i], cs_A->x[i]);
                i++;
                cs_A->i[i] = negative - 1;
                cs_A->p[i] = data + Vcnt;
                cs_A->x[i] = -1;
                //printf("[%d, %d]: %.4f\n", cs_A->i[i], cs_A->p[i], cs_A->x[i]);
                i++;
            } else if (!negative) {
                cs_A->i[i] = data + Vcnt;
                cs_A->p[i] = positive - 1;
                cs_A->x[i] = 1;
                //printf("[%d, %d]: %.4f\n", cs_A->i[i], cs_A->p[i], cs_A->x[i]);
                i++;
                cs_A->i[i] = positive - 1;
                cs_A->p[i] = data + Vcnt;
                cs_A->x[i] = 1;
                //printf("[%d, %d]: %.4f\n", cs_A->i[i], cs_A->p[i], cs_A->x[i]);
                i++;
            } else {
                cs_A->i[i] = data + Vcnt;
                cs_A->p[i] = positive - 1;
                cs_A->x[i] = 1;
                //printf("[%d, %d]: %.4f\n", cs_A->i[i], cs_A->p[i], cs_A->x[i]);
                i++;
                cs_A->i[i] = data + Vcnt;
                cs_A->p[i] = negative - 1;
                cs_A->x[i] = -1;
                //printf("[%d, %d]: %.4f\n", cs_A->i[i], cs_A->p[i], cs_A->x[i]);
                i++;
                cs_A->i[i] = positive - 1;
                cs_A->p[i] = data + Vcnt;
                cs_A->x[i] = 1;
                //printf("[%d, %d]: %.4f\n", cs_A->i[i], cs_A->p[i], cs_A->x[i]);
                i++;
                cs_A->i[i] = negative - 1;
                cs_A->p[i] = data + Vcnt;
                cs_A->x[i] = -1;
                //printf("[%d, %d]: %.4f\n", cs_A->i[i], cs_A->p[i], cs_A->x[i]);
                i++;
            }
            cs_b[data + Vcnt] = atof(n->value);
            Vcnt++;
        }
        n = n->nxt;
    }
    if (todoCSCG || todoCSBiCG) {
        for (i = data; i < data + a2; i++) {
            gsl_matrix_set(diag_A, i, i, 1);
        }
    }
    cs_C = cs_compress(cs_A);
    if (cs_C == NULL) printf("compress fail\n");
    else printf("compress success\n");
    if (cs_dupl(cs_C)) printf("dupl success\n");
}

void populateTransientSystem(circuitList *n, HASHTBL *hashtbl, gsl_matrix * tilda_G, gsl_matrix * tilda_C, gsl_vector * b) {
    //printf("populate trans\n");
    int positive;
    int negative;
    float prevVal;
    int Vcnt = 0;
    int Lcnt = 0;
    while (n != NULL) {
        //printf("trans1\n");
        positive = atoi(hashtbl_get(hashtbl, n->positive));
        negative = atoi(hashtbl_get(hashtbl, n->negative));
        if (strstr(n->name, "R") || strstr(n->name, "r")) {
            if (!positive) {
                prevVal = gsl_matrix_get(tilda_G, negative - 1, negative - 1);
                gsl_matrix_set(tilda_G, negative - 1, negative - 1, 1 / atof(n->value) + prevVal);
            } else if (!negative) {
                prevVal = gsl_matrix_get(tilda_G, positive - 1, positive - 1);
                gsl_matrix_set(tilda_G, positive - 1, positive - 1, 1 / atof(n->value) + prevVal);
            } else {
                prevVal = gsl_matrix_get(tilda_G, positive - 1, positive - 1);
                gsl_matrix_set(tilda_G, positive - 1, positive - 1, 1 / atof(n->value) + prevVal);
                prevVal = gsl_matrix_get(tilda_G, positive - 1, negative - 1);
                gsl_matrix_set(tilda_G, positive - 1, negative - 1, -(1 / atof(n->value)) + prevVal);
                prevVal = gsl_matrix_get(tilda_G, negative - 1, positive - 1);
                gsl_matrix_set(tilda_G, negative - 1, positive - 1, -(1 / atof(n->value)) + prevVal);
                prevVal = gsl_matrix_get(tilda_G, negative - 1, negative - 1);
                gsl_matrix_set(tilda_G, negative - 1, negative - 1, 1 / atof(n->value) + prevVal);
            }
        } else if (strstr(n->name, "I") || strstr(n->name, "i")) {
            if (!positive) {
                prevVal = gsl_vector_get(tran_b, negative - 1);
                gsl_vector_set(tran_b, negative - 1, atof(n->value) + prevVal);
            } else if (!negative) {
                prevVal = gsl_vector_get(tran_b, positive - 1);
                gsl_vector_set(tran_b, positive - 1, -atof(n->value) + prevVal);
            } else {
                prevVal = gsl_vector_get(tran_b, positive - 1);
                gsl_vector_set(tran_b, positive - 1, -atof(n->value) + prevVal);
                prevVal = gsl_vector_get(tran_b, negative - 1);
                gsl_vector_set(tran_b, negative - 1, atof(n->value) + prevVal);
            }
        } else if (strstr(n->name, "V") || strstr(n->name, "v")) {
            if (!positive) {
                prevVal = gsl_matrix_get(tilda_G, data + Vcnt, negative - 1);
                gsl_matrix_set(tilda_G, data + Vcnt, negative - 1, -1 + prevVal);
                prevVal = gsl_matrix_get(tilda_G, negative - 1, data + Vcnt);
                gsl_matrix_set(tilda_G, negative - 1, data + Vcnt, -1 + prevVal);
                gsl_vector_set(tran_b, data + Vcnt, atof(n->value));
            } else if (!negative) {
                prevVal = gsl_matrix_get(tilda_G, data + Vcnt, positive - 1);
                gsl_matrix_set(tilda_G, data + Vcnt, positive - 1, 1 + prevVal);
                prevVal = gsl_matrix_get(tilda_G, positive - 1, data + Vcnt);
                gsl_matrix_set(tilda_G, positive - 1, data + Vcnt, 1 + prevVal);
                gsl_vector_set(tran_b, data + Vcnt, atof(n->value));
            } else {
                prevVal = gsl_matrix_get(tilda_G, data + Vcnt, positive - 1);
                gsl_matrix_set(tilda_G, data + Vcnt, positive - 1, 1 + prevVal);
                prevVal = gsl_matrix_get(tilda_G, data + Vcnt, negative - 1);
                gsl_matrix_set(tilda_G, data + Vcnt, negative - 1, -1 + prevVal);
                prevVal = gsl_matrix_get(tilda_G, positive - 1, data + Vcnt);
                gsl_matrix_set(tilda_G, positive - 1, data + Vcnt, 1 + prevVal);
                prevVal = gsl_matrix_get(tilda_G, negative - 1, data + Vcnt);
                gsl_matrix_set(tilda_G, negative - 1, data + Vcnt, -1 + prevVal);
                gsl_vector_set(tran_b, data + Vcnt, atof(n->value));
            }
            Vcnt++;
        } else if (strstr(n->name, "L") || strstr(n->name, "l")) {
            if (!positive) {
                prevVal = gsl_matrix_get(tilda_G, data + Lcnt, negative - 1);
                gsl_matrix_set(tilda_G, data + Lcnt, negative - 1, -1 + prevVal);
                prevVal = gsl_matrix_get(tilda_G, negative - 1, data + Lcnt);
                gsl_matrix_set(tilda_G, negative - 1, data + Lcnt, -1 + prevVal);


            } else if (!negative) {
                prevVal = gsl_matrix_get(tilda_G, data + Lcnt, positive - 1);
                gsl_matrix_set(tilda_G, data + Lcnt, positive - 1, 1 + prevVal);
                prevVal = gsl_matrix_get(tilda_G, positive - 1, data + Lcnt);
                gsl_matrix_set(tilda_G, positive - 1, data + Lcnt, 1 + prevVal);

            } else {
                prevVal = gsl_matrix_get(tilda_G, data + Lcnt, positive - 1);
                gsl_matrix_set(tilda_G, data + Lcnt, positive - 1, 1 + prevVal);
                prevVal = gsl_matrix_get(tilda_G, data + Lcnt, negative - 1);
                gsl_matrix_set(tilda_G, data + Lcnt, negative - 1, -1 + prevVal);
                prevVal = gsl_matrix_get(tilda_G, positive - 1, data + Lcnt);
                gsl_matrix_set(tilda_G, positive - 1, data + Lcnt, 1 + prevVal);
                prevVal = gsl_matrix_get(tilda_G, negative - 1, data + Lcnt);
                gsl_matrix_set(tilda_G, negative - 1, data + Lcnt, -1 + prevVal);

            }
            gsl_matrix_set(tilda_C, data + Lcnt, data + Lcnt, -atof(n->value));
            Lcnt++;

        } else if (strstr(n->name, "C") || strstr(n->name, "c")) {
            if (!positive) {
                prevVal = gsl_matrix_get(tilda_C, negative - 1, negative - 1);
                gsl_matrix_set(tilda_C, negative - 1, negative - 1, atof(n->value) + prevVal);
            } else if (!negative) {
                prevVal = gsl_matrix_get(tilda_C, positive - 1, positive - 1);
                gsl_matrix_set(tilda_C, positive - 1, positive - 1, atof(n->value) + prevVal);
            } else {
                prevVal = gsl_matrix_get(tilda_C, positive - 1, positive - 1);
                gsl_matrix_set(tilda_C, positive - 1, positive - 1, atof(n->value) + prevVal);
                prevVal = gsl_matrix_get(tilda_C, positive - 1, negative - 1);
                gsl_matrix_set(tilda_C, positive - 1, negative - 1, -atof(n->value) + prevVal);
                prevVal = gsl_matrix_get(tilda_C, negative - 1, positive - 1);
                gsl_matrix_set(tilda_C, negative - 1, positive - 1, -atof(n->value) + prevVal);
                prevVal = gsl_matrix_get(tilda_C, negative - 1, negative - 1);
                gsl_matrix_set(tilda_C, negative - 1, negative - 1, atof(n->value) + prevVal);
            }
        }
        n = n->nxt;
    }
    //printf("exit\n");
    //printGSLMatrix("tilda_G", tilda_G);
    //printf("\n");
    //printGSLMatrix("tilda_C", tilda_C);
    //printf("\n");
    //printGslVector("tran_b", tran_b);
}

void calcNonZeroElem(circuitList *n, HASHTBL *hashtbl) {
    int positive;
    int negative;
    nonZeroElem = 0;
    while (n != NULL) {
        positive = atoi(hashtbl_get(hashtbl, n->positive));
        negative = atoi(hashtbl_get(hashtbl, n->negative));
        if (strstr(n->name, "R") || strstr(n->name, "r")) {
            if (!positive || !negative) nonZeroElem++;
            else nonZeroElem += 4;
        } else if (strstr(n->name, "V") || strstr(n->name, "v")) {
            if (!positive || !negative) nonZeroElem += 2;
            else nonZeroElem += 4;
        }
        n = n->nxt;
    }
}

void populateB(circuitList *n, gsl_vector *tran, float i, int k, float prevVal, float value) {
    int positive, negative;
    int Vcnt = 0;
    while (n != NULL) {
        //printf("transient_spec=%s\n", n->transient_spec);
        positive = atoi(hashtbl_get(hashtbl, n->positive));
        //printf("positive=%d\n", positive);
        negative = atoi(hashtbl_get(hashtbl, n->negative));
        //printf("negative=%d\n", negative);
        if (strstr(n->name, "I") || strstr(n->name, "i")) {
            //printf("if_i\n");
            if (!positive && strstr(n->transient_spec, "EXP")) {
                //printf("EXP\n");
                prevVal = gsl_vector_get(tran, negative - 1);
                value = exp_solve(i, n->exp);
                gsl_vector_set(tran, negative - 1, value + prevVal);
            } else if (!positive && strstr(n->transient_spec, "SIN")) {
                //printf("SIN\n");
                prevVal = gsl_vector_get(tran, negative - 1);
                value = sin_solve(i, n->sin);
                gsl_vector_set(tran, negative - 1, value + prevVal);
            } else if (!positive && strstr(n->transient_spec, "PULSE")) {
                //printf("PULSE\n");
                prevVal = gsl_vector_get(tran, negative - 1);
                value = pulse_solve(k, i, n->pulse);
                gsl_vector_set(tran, negative - 1, value + prevVal);
            } else if (!positive && strstr(n->transient_spec, "PWL")) {
                //n->pwl = pwlRoot;
                //printf("PWL1\n");
                prevVal = gsl_vector_get(tran, negative - 1);
                value = pwl_solve(i, pwlRoot);
                //printf("value=%f\n", value);
                gsl_vector_set(tran, negative - 1, value + prevVal);
            } else if (!negative && strstr(n->transient_spec, "SIN")) {
                //printf("SIN\n");
                prevVal = gsl_vector_get(tran, positive - 1);
                value = sin_solve(i, n->sin);
                gsl_vector_set(tran, positive - 1, -value + prevVal);
            } else if (!negative && strstr(n->transient_spec, "PULSE")) {
                //printf("PULSE\n");
                prevVal = gsl_vector_get(tran, positive - 1);
                value = pulse_solve(k, i, n->pulse);
                gsl_vector_set(tran, positive - 1, -value + prevVal);
            } else if (!negative && strstr(n->transient_spec, "PWL")) {
                //n->pwl = pwlRoot;
                //printf("PWL2\n");
                prevVal = gsl_vector_get(tran, positive - 1);
                value = pwl_solve(i, pwlRoot);
                gsl_vector_set(tran, positive - 1, -value + prevVal);
            } else if (!negative && strstr(n->transient_spec, "EXP")) {
                //printf("EXP\n");
                prevVal = gsl_vector_get(tran, positive - 1);
                value = exp_solve(i, n->exp);
                gsl_vector_set(tran, positive - 1, -value + prevVal);
            } else if (strstr(n->transient_spec, "EXP")) {
                //printf("EXP\n");
                prevVal = gsl_vector_get(tran, positive - 1);
                value = exp_solve(i, n->exp);
                gsl_vector_set(tran, positive - 1, -value + prevVal);
                prevVal = gsl_vector_get(tran, negative - 1);
                value = exp_solve(i, n->exp);
                gsl_vector_set(tran, negative - 1, value + prevVal);
            } else if (strstr(n->transient_spec, "SIN")) {
                //printf("SIN\n");
                prevVal = gsl_vector_get(tran, positive - 1);
                value = sin_solve(i, n->sin);
                gsl_vector_set(tran, positive - 1, -value + prevVal);
                prevVal = gsl_vector_get(tran, negative - 1);
                value = sin_solve(i, n->sin);
                gsl_vector_set(tran, negative - 1, value + prevVal);
            } else if (strstr(n->transient_spec, "PULSE")) {
                //printf("PULSE\n");
                prevVal = gsl_vector_get(tran, positive - 1);
                value = pulse_solve(k, i, n->pulse);
                gsl_vector_set(tran, positive - 1, -value + prevVal);
                prevVal = gsl_vector_get(tran, negative - 1);
                value = pulse_solve(k, i, n->pulse);
                gsl_vector_set(tran, negative - 1, value + prevVal);
            } else if (strstr(n->transient_spec, "PWL")) {
                n->pwl = pwlRoot;
                //printf("PWL3\n");
                prevVal = gsl_vector_get(tran, positive - 1);
                value = pwl_solve(i, n->pwl);
                gsl_vector_set(tran, positive - 1, -value + prevVal);
                prevVal = gsl_vector_get(tran, negative - 1);
                n->pwl = pwlRoot;
                value = pwl_solve(i, n->pwl);
                gsl_vector_set(tran, negative - 1, value + prevVal);
            } else if (!positive && (strstr(n->transient_spec, "NO"))) {
                //printf("NO\n");
                prevVal = gsl_vector_get(tran, negative - 1);
                gsl_vector_set(tran, negative - 1, atof(n->value) + prevVal);
            } else if (!negative && (strstr(n->transient_spec, "NO"))) {
                //printf("NO\n");
                prevVal = gsl_vector_get(tran, positive - 1);
                gsl_vector_set(tran, positive - 1, -atof(n->value) + prevVal);
            } else if (strstr(n->transient_spec, "NO")) {
                //printf("NO\n");
                prevVal = gsl_vector_get(tran, positive - 1);
                gsl_vector_set(tran, positive - 1, -atof(n->value) + prevVal);
                prevVal = gsl_vector_get(tran, negative - 1);
                gsl_vector_set(tran, negative - 1, atof(n->value) + prevVal);
            }
        } else if (strstr(n->name, "V") || strstr(n->name, "v")) {
            //printf("if_v\n");
            //printf("transient_spec=%s\n", n->transient_spec);
            if (strstr(n->transient_spec, "EXP")) {
                //printf("EXP\n");
                value = exp_solve(i, n->exp);
                //printf("value=%f\n", value);
                gsl_vector_set(tran, data + Vcnt, value);
            } else if (strstr(n->transient_spec, "SIN")) {
                //printf("SIN\n");
                value = sin_solve(i, n->sin);
                //printf("value=%f\n", value);
                gsl_vector_set(tran, data + Vcnt, value);
            } else if (strstr(n->transient_spec, "PULSE")) {
                //printf("PULSE\n");
                value = pulse_solve(k, i, n->pulse);
                //printf("value=%f\n", value);
                gsl_vector_set(tran, data + Vcnt, value);
            } else if (strstr(n->transient_spec, "PWL")) {
                n->pwl = pwlRoot;
                //printf("PWL1\n");
                value = pwl_solve(i, n->pwl);
                //printf("value=%f\n", value);
                gsl_vector_set(tran, data + Vcnt, value);
            } else if (strstr(n->transient_spec, "NO")) {
                //printf("NO\n");
                gsl_vector_set(tran, data + Vcnt, atof(n->value));
            }
            Vcnt++;
        } else {
        }
        n = n->nxt;
    }
}

void solveSystem() {
    if (todoCholesky) {
        printf("\n*** Cholesky ***\n");
        gsl_matrix *LLU = gsl_matrix_calloc(data, data);
        gsl_matrix_memcpy(LLU, A);
        gsl_vector *b_temp = gsl_vector_calloc(data);
        gsl_vector_memcpy(b_temp, b);
        gsl_linalg_cholesky_decomp(LLU);
        if (hasDC) {
            circuitList *node = circuitListGet(dc->value);
            int positive = atoi(hashtbl_get(hashtbl, node->positive));
            int negative = atoi(hashtbl_get(hashtbl, node->negative));
            float i;
            for (i = atof(dc->start); i <= atof(dc->end); i += atof(dc->step)) {
                
                if (positive != 0) gsl_vector_set(b, positive - 1, -i);
                if (negative != 0) gsl_vector_set(b, negative - 1, i);
                //printGslVector("b", b);
                gsl_linalg_cholesky_solve(LLU, b, x);
                if (todoPlot) printGslPlot(plotRoot, x);
                else printGslVector("x", x);
            }
        }
        gsl_linalg_cholesky_solve(LLU, b_temp, x);
        if (todoPlot) printGslPlot(plotRoot, x);
        else printGslVector("x", x);
    } else if (todoCG) {
        printf("\n*** CG ***\n");
        gsl_matrix *cg_A = gsl_matrix_calloc(data + a2, data + a2);
        gsl_matrix_memcpy(cg_A, A);
        gsl_matrix *cg_M = gsl_matrix_calloc(data + a2, data + a2);
        gsl_vector *cg_p = gsl_vector_calloc(data + a2);
        gsl_vector *cg_q = gsl_vector_calloc(data + a2);
        gsl_vector *cg_r = gsl_vector_calloc(data + a2);
        gsl_vector *cg_x = gsl_vector_calloc(data + a2);
        gsl_vector *cg_z = gsl_vector_calloc(data + a2);
        double cg_alpha;
        double cg_beta;
        double cg_rh0;
        double cg_rh1;
        gsl_vector_memcpy(cg_r, b);
        int iter = 0;
        cg_M = composeGslDiagonal(A);
        while (gsl_blas_dnrm2(cg_r) / gsl_blas_dnrm2(b) > itol && iter < data + a2) {
            iter++;
            cg_z = solveGslPreconditioner(cg_M, cg_z, cg_r);
            gsl_blas_ddot(cg_r, cg_z, &cg_rh0);
            if (iter == 1) gsl_vector_memcpy(cg_p, cg_z);
            else {
                cg_beta = cg_rh0 / cg_rh1;
                gsl_vector_scale(cg_p, cg_beta);
                gsl_vector_add(cg_p, cg_z);
            }
            cg_rh1 = cg_rh0;
            gsl_blas_dgemv(CblasNoTrans, 1.0, cg_A, cg_p, 0.0, cg_q);
            gsl_blas_ddot(cg_p, cg_q, &cg_alpha);
            cg_alpha = cg_rh0 / cg_alpha;
            gsl_blas_daxpy(cg_alpha, cg_p, cg_x);
            gsl_blas_daxpy(-cg_alpha, cg_q, cg_r);
        }

        if (todoPlot) printGslPlot(plotRoot, cg_x);
        else printGslVector("x", cg_x);
    } else if (todoBiCG) {
        printf("\n*** BiCG ***\n");
        gsl_matrix *cg_A = gsl_matrix_calloc(data + a2, data + a2);
        gsl_matrix_memcpy(cg_A, A);
        gsl_matrix *bicg_M = gsl_matrix_calloc(data + a2, data + a2);
        gsl_vector *bicg_p = gsl_vector_calloc(data + a2);
        gsl_vector *bicg_pt = gsl_vector_calloc(data + a2);
        gsl_vector *bicg_q = gsl_vector_calloc(data + a2);
        gsl_vector *bicg_qt = gsl_vector_calloc(data + a2);
        gsl_vector *bicg_r = gsl_vector_calloc(data + a2);
        gsl_vector *bicg_rt = gsl_vector_calloc(data + a2);
        gsl_vector *bicg_x = gsl_vector_calloc(data + a2);
        gsl_vector *bicg_z = gsl_vector_calloc(data + a2);
        gsl_vector *bicg_zt = gsl_vector_calloc(data + a2);
        double bicg_alpha;
        double bicg_beta;
        double omega;
        double bicg_rh0;
        double bicg_rh1;
        int iter = 0;
        gsl_vector_memcpy(bicg_r, b);
        gsl_vector_memcpy(bicg_rt, b);
        bicg_M = composeGslDiagonal(A);
        while (gsl_blas_dnrm2(bicg_r) / gsl_blas_dnrm2(b) > itol && iter < data + a2) {
            iter++;
            bicg_z = solveGslPreconditioner(bicg_M, bicg_z, bicg_r);
            bicg_zt = solveGslPreconditioner(bicg_M, bicg_zt, bicg_rt);
            gsl_blas_ddot(bicg_z, bicg_rt, &bicg_rh0);
            if (fabs(bicg_rh0) < atof("1e-14")) {
                printf("algorithm failure: rh0\n");
                break;
            }
            if (iter == 1) {
                gsl_vector_memcpy(bicg_p, bicg_z);
                gsl_vector_memcpy(bicg_pt, bicg_zt);
            } else {
                bicg_beta = bicg_rh0 / bicg_rh1;
                gsl_vector_scale(bicg_p, bicg_beta);
                gsl_vector_add(bicg_p, bicg_z);
                gsl_vector_scale(bicg_pt, bicg_beta);
                gsl_vector_add(bicg_pt, bicg_zt);
            }
            bicg_rh1 = bicg_rh0;
            gsl_blas_dgemv(CblasNoTrans, 1.0, cg_A, bicg_p, 0.0, bicg_q);
            gsl_blas_dgemv(CblasTrans, 1.0, cg_A, bicg_pt, 0.0, bicg_qt);
            gsl_blas_ddot(bicg_pt, bicg_q, &omega);
            if (fabs(omega) < atof("1e-14")) {
                printf("algorithm failure: omega\n");
                break;
            }
            bicg_alpha = bicg_rh0 / omega;
            gsl_blas_daxpy(bicg_alpha, bicg_p, bicg_x);
            gsl_blas_daxpy(-bicg_alpha, bicg_q, bicg_r);
            gsl_blas_daxpy(-bicg_alpha, bicg_qt, bicg_rt);
        }
        if (todoPlot) printGslPlot(plotRoot, bicg_x);
        else printGslVector("x", bicg_x);
    } else {
        printf("\n*** LU ***\n");
        gsl_permutation *p = gsl_permutation_alloc(data + a2);
        gsl_matrix *PA = gsl_matrix_calloc(data + a2, data + a2);
        gsl_matrix_memcpy(PA, A);
        int s;
        gsl_linalg_LU_decomp(PA, p, &s);
        if (hasDC) {
            if (dc->value[0] == 'v' || dc->value[0] == 'V') {
                char *dcValue = (char *) malloc(sizeof (dc->value));
                strcpy(dcValue, dc->value);
                char *token = strtok(dcValue, "Vv");
                int position = atoi(token);
                float i;
                for (i = atof(dc->start); i <= atof(dc->end); i += atof(dc->step)) {
                    gsl_vector_set(b, data + position - 1, i);
                    //printGslVector("b", b);
                    gsl_linalg_LU_solve(PA, p, b, x);

                    if (todoPlot) printGslPlot(plotRoot, x);
                    else printGslVector("x", x);
                }
            } else {
                circuitList *node = circuitListGet(dc->value);
                int positive = atoi(hashtbl_get(hashtbl, node->positive));
                int negative = atoi(hashtbl_get(hashtbl, node->negative));
                float i;
                for (i = atof(dc->start); i <= atof(dc->end); i += atof(dc->step)) {
                    if (positive != 0) gsl_vector_set(b, positive - 1, -i);
                    if (negative != 0) gsl_vector_set(b, negative - 1, i);
                    //printGslVector("b", b);
                    gsl_linalg_LU_solve(PA, p, b, x);
                    if (todoPlot) printGslPlot(plotRoot, x);
                    else printGslVector("x", x);

                }
            }
        } else {
            gsl_linalg_LU_solve(PA, p, b, x);

            if (todoPlot) printGslPlot(plotRoot, x);
            else printGslVector("x", x);
        }
    }
}

void solveCSSystem() {
    if (todoCSLU) {
        printf("\n*** CSsparse LU ***\n");
        css *cs_S = cs_sqr(2, cs_C, 0);
        csn *cs_N = cs_lu(cs_C, cs_S, 1);
        if (cs_N == NULL) printf("lu fail\n");
        else printf("lu success\n");
        if (cs_ipvec(cs_N->pinv, cs_b, cs_x, data + a2)) printf("ipvec success\n");
        else printf("ipvec failure\n");
        if (cs_lsolve(cs_N->L, cs_x)) printf("lsolve success\n");
        else printf("lsolve failure\n");
        if (cs_usolve(cs_N->U, cs_x)) printf("usolve success\n");
        else printf("usolve failure\n");
        if (cs_ipvec(cs_S->q, cs_x, cs_b, data + a2)) printf("ipvec success\n");
        else printf("ipvec failure\n");
        if (todoPlot) printPlot(plotRoot, cs_b);
        else printVector("x", cs_b, data + a2);
    } else if (todoCSCholesky) {
        printf("\n*** CSsparse Cholesky ***\n");
        css *cs_S = cs_schol(1, cs_C);
        csn *cs_N = cs_chol(cs_C, cs_S);
        if (cs_S == NULL) printf("schol fail\n");
        else printf("schol success\n");
        if (cs_N == NULL) printf("chol fail\n");
        else printf("chol success\n");
        if (cs_ipvec(cs_S->pinv, cs_b, cs_x, data + a2)) printf("ipvec success\n");
        else printf("ipvec failure\n");
        if (cs_lsolve(cs_N->L, cs_x)) printf("lsolve success\n");
        else printf("lsolve failure\n");
        if (cs_ltsolve(cs_N->L, cs_x)) printf("ltsolve success\n");
        else printf("ltsolve failure\n");
        if (cs_pvec(cs_S->pinv, cs_x, cs_b, data + a2)) printf("pvec success\n");
        else printf("pvec failure\n");
        if (todoPlot) printPlot(plotRoot, cs_b);
        else printVector("x", cs_b, data + a2);
    } else if (todoCSCG) {
        printf("\n*** CSsparse CG ***\n");
        gsl_matrix *cg_M = gsl_matrix_calloc(data + a2, data + a2);
        gsl_vector *cg_b = doubleToGslVector(cs_b, data + a2);
        gsl_vector *cg_p = gsl_vector_calloc(data + a2);
        gsl_vector *cg_r = gsl_vector_calloc(data + a2);
        gsl_vector *cg_x = gsl_vector_calloc(data + a2);
        gsl_vector *cg_z = gsl_vector_calloc(data + a2);
        double *cg_q;
        double cg_alpha;
        double cg_beta;
        double cg_rh0;
        double cg_rh1;
        int iter = 0;
        gsl_vector_memcpy(cg_r, cg_b);
        gsl_matrix_memcpy(cg_M, diag_A);
        while (gsl_blas_dnrm2(cg_r) / gsl_blas_dnrm2(cg_b) > itol && iter < data + a2) {
            iter++;
            cg_z = solveGslPreconditioner(cg_M, cg_z, cg_r);
            gsl_blas_ddot(cg_r, cg_z, &cg_rh0);
            if (iter == 1) gsl_vector_memcpy(cg_p, cg_z);
            else {
                cg_beta = cg_rh0 / cg_rh1;
                gsl_vector_scale(cg_p, cg_beta);
                gsl_vector_add(cg_p, cg_z);
            }
            cg_rh1 = cg_rh0;
            cg_q = cs_gaxp(cs_C, gslVectorToDouble(cg_p));
            gsl_blas_ddot(cg_p, doubleToGslVector(cg_q, data + a2), &cg_alpha);
            cg_alpha = cg_rh0 / cg_alpha;
            gsl_blas_daxpy(cg_alpha, cg_p, cg_x);
            gsl_blas_daxpy(-cg_alpha, doubleToGslVector(cg_q, data + a2), cg_r);
        }
        if (todoPlot) printGslPlot(plotRoot, cg_x);
        else printGslVector("x", cg_x);
    } else if (todoCSBiCG) {
        printf("\n*** CSsparse BiCG ***\n");
        gsl_matrix *bicg_M = gsl_matrix_calloc(data + a2, data + a2);
        gsl_vector *bicg_b = doubleToGslVector(cs_b, data + a2);
        gsl_vector *bicg_p = gsl_vector_calloc(data + a2);
        gsl_vector *bicg_pt = gsl_vector_calloc(data + a2);
        gsl_vector *bicg_r = gsl_vector_calloc(data + a2);
        gsl_vector *bicg_rt = gsl_vector_calloc(data + a2);
        gsl_vector *bicg_x = gsl_vector_calloc(data + a2);
        gsl_vector *bicg_z = gsl_vector_calloc(data + a2);
        gsl_vector *bicg_zt = gsl_vector_calloc(data + a2);
        double *bicg_q;
        double *bicg_qt;
        double bicg_alpha;
        double bicg_beta;
        double omega;
        double bicg_rh0;
        double bicg_rh1;
        int iter = 0;
        gsl_vector_memcpy(bicg_r, bicg_b);
        gsl_vector_memcpy(bicg_rt, bicg_b);
        cs *cs_Ct = cs_transpose(cs_C, 1);
        gsl_matrix_memcpy(bicg_M, diag_A);
        while (gsl_blas_dnrm2(bicg_r) / gsl_blas_dnrm2(bicg_b) > itol && iter < 1000) {
            iter++;
            bicg_z = solveGslPreconditioner(bicg_M, bicg_z, bicg_r);
            bicg_zt = solveGslPreconditioner(bicg_M, bicg_zt, bicg_rt);
            gsl_blas_ddot(bicg_z, bicg_rt, &bicg_rh0);
            if (fabs(bicg_rh0) < atof("1e-14")) {
                printf("algorithm failure: rh0\n");
                break;
            }
            if (iter == 1) {
                gsl_vector_memcpy(bicg_p, bicg_z);
                gsl_vector_memcpy(bicg_pt, bicg_zt);
            } else {
                bicg_beta = bicg_rh0 / bicg_rh1;
                gsl_vector_scale(bicg_p, bicg_beta);
                gsl_vector_add(bicg_p, bicg_z);
                gsl_vector_scale(bicg_pt, bicg_beta);
                gsl_vector_add(bicg_pt, bicg_zt);
            }
            bicg_rh1 = bicg_rh0;
            bicg_q = cs_gaxp(cs_C, gslVectorToDouble(bicg_p));
            bicg_qt = cs_gaxp(cs_Ct, gslVectorToDouble(bicg_pt));
            gsl_blas_ddot(bicg_pt, doubleToGslVector(bicg_q, data + a2), &omega);
            if (fabs(omega) < atof("1e-14")) {
                printf("algorithm failure: omega\n");
                break;
            }
            bicg_alpha = bicg_rh0 / omega;
            gsl_blas_daxpy(bicg_alpha, bicg_p, bicg_x);
            gsl_blas_daxpy(-bicg_alpha, doubleToGslVector(bicg_q, data + a2), bicg_r);
            gsl_blas_daxpy(-bicg_alpha, doubleToGslVector(bicg_qt, data + a2), bicg_rt);
        }
        //printf("iter: %d\n", iter);
        if (todoPlot) printGslPlot(plotRoot, bicg_x);
        else printGslVector("x", bicg_x);
    }
}

void solveTransientSystem(circuitList *n) {
    float i;
    int k = 0;
    int positive;
    int negative;
    int Vcnt = 0;
    float prevVal;
    float value;
    float h = atof(tran->step);
    //printf("h=%f\n", h);
    if (todoBE) {
        printf("\n*** Backward Euler***\n");
        gsl_matrix *tilda_G_be = gsl_matrix_alloc(data + a2, data + a2);
        gsl_matrix_memcpy(tilda_G_be, tilda_G);
        gsl_matrix *tilda_C_be = gsl_matrix_alloc(data + a2, data + a2);
        gsl_matrix_memcpy(tilda_C_be, tilda_C);
        gsl_matrix_scale(tilda_C_be, 1 / h);
        gsl_matrix_add(tilda_G_be, tilda_C_be);
        LU_solve(tilda_G, tran_x, tran_b);
        gsl_blas_dgemv(CblasNoTrans, 1.0, tilda_C_be, tran_x, 0.0, tran_x);
        for (i = h; i <= atof(tran->end); i += h) {
            //printf("i=%f\n", i);
            k++;
            populateB(n, tran_b, i, k, prevVal, value);
            //printf("exit\n");
            gsl_vector_add(tran_b, tran_x);
            LU_solve(tilda_G_be, tran_x, tran_b);

            //printf("\n");
            if (todoPlot) printGslPlot(plotRoot, tran_x);
            else printGslVector("tran_x", tran_x);
            //printf("\n");
            n = circuitRoot;
            gsl_blas_dgemv(CblasNoTrans, 1.0, tilda_C_be, tran_x, 0.0, tran_x);
        }
    } else if (todoTransient) {
        printf("\n*** Trapezodial***\n");
        gsl_matrix *tilda_G_tr = gsl_matrix_alloc(data + a2, data + a2);
        gsl_matrix_memcpy(tilda_G_tr, tilda_G);
        gsl_matrix *tilda_C_tr = gsl_matrix_alloc(data + a2, data + a2);
        gsl_matrix *tilda_C_2h = gsl_matrix_alloc(data + a2, data + a2);
        gsl_matrix_memcpy(tilda_C_2h, tilda_C);
        gsl_matrix_scale(tilda_C_2h, 2 / h);
        gsl_matrix_memcpy(tilda_C_tr, tilda_C_2h);
        gsl_matrix_add(tilda_G_tr, tilda_C_2h);
        LU_solve(tilda_G, tran_x, tran_b);
        gsl_matrix_add(tilda_C_2h, tilda_G);
        gsl_matrix_sub(tilda_G_tr, tilda_C_tr);
        gsl_blas_dgemv(CblasNoTrans, 1.0, tilda_G_tr, tran_x, 0.0, tran_x);
        //printf("transient_spec=%s\n", n->transient_spec);
        positive = atoi(hashtbl_get(hashtbl, n->positive));
        //printf("positive=%d\n", positive);
        negative = atoi(hashtbl_get(hashtbl, n->negative));
        //printf("negative=%d\n", negative);
        populateB(n, tran_b, i, k, prevVal, value);
        n = circuitRoot;
        k++;
        for (i = h; i <= atof(tran->end); i += h) {
            Vcnt = 0;
            populateB(n, tran_b_nxt, i, k, prevVal, value);
            gsl_vector_add(tran_b, tran_b_nxt);
            gsl_vector_sub(tran_b, tran_x);
            LU_solve(tilda_C_2h, tran_x, tran_b);
            k++;
            //printf("\n");
            if (todoPlot) printGslPlot(plotRoot, tran_x);
            else printGslVector("tran_x", tran_x);
            //printf("\n");
            n = circuitRoot;
            gsl_vector_memcpy(tran_b, tran_b_nxt);
            gsl_blas_dgemv(CblasNoTrans, 1.0, tilda_G_tr, tran_x, 0.0, tran_x);
        }
    }
}