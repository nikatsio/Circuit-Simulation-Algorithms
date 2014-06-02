#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdbool.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include "csparse.h"
#include "list.h"
#include "hashtbl.h"
#include "parser.h"
#include "table.h"
#include "util.h"

circuitList *circuitRoot;
circuitList *circuitCurr;
mosList *mosRoot;
mosList *mosCurr;
bjtList *bjtRoot;
bjtList *bjtCurr;
dcList *dc;
tranList *tran;
plotList *plotRoot;
plotList *plotCurr;
pwlList *pwlRoot;
pwlList *pwlCurr;

int element = 0;
bool hasDC = false;

void listInit(char *list) {
    //printf("Initializing %s\n", list);
    if (strcmp(list, "circuit") == 0) {
        circuitRoot = NULL;
        circuitCurr = NULL;
    } else if (strcmp(list, "mos") == 0) {
        mosRoot = NULL;
        mosCurr = NULL;
    } else if (strcmp(list, "bjt") == 0) {
        bjtRoot = NULL;
        bjtCurr = NULL;
    } else if (strcmp(list, "plot") == 0) {
        plotRoot = NULL;
        plotCurr = NULL;
    }
}

int listHasElement(char *list) {
    if (strcmp(list, "circuit") == 0) return (circuitRoot != NULL);
    else if (strcmp(list, "mos") == 0) return (mosRoot != NULL);
    else if (strcmp(list, "bjt") == 0) return (bjtRoot != NULL);
    else if (strcmp(list, "plot") == 0) return (plotRoot != NULL);
    else if (strcmp(list, "pwl") == 0) return (pwlRoot != NULL);
}

void circuitListAdd(circuitList *node) {
    if (!listHasElement("circuit")) {
        circuitRoot = node;
        circuitCurr = node;
    } else {
        circuitCurr->nxt = node;
        circuitCurr = node;
    }
    element++;
}

void *circuitListGet(char *name) {
    //circuitList *node = (circuitList *) malloc(sizeof (circuitList));
    circuitList *node = NULL;
    node = circuitRoot;
    if (listHasElement("circuit")) {
        while (node) {
            if (!strcmp(node->name, name)) return node;
            node = node->nxt;
        }
    } else return node;
}

void mosListAdd(mosList *node) {
    if (!listHasElement("mos")) {
        mosRoot = node;
        mosCurr = node;
    } else {
        mosCurr->nxt = node;
        mosCurr = node;
    }
}

void bjtListAdd(bjtList *node) {
    if (!listHasElement("bjt")) {
        bjtRoot = node;
        bjtCurr = node;
    } else {
        bjtCurr->nxt = node;
        bjtCurr = node;
    }
}

void plotListAdd(plotList *node) {
    //printf("plotListAdd: %s\n", node->value);
    if (!listHasElement("plot")) {
        plotRoot = node;
        plotCurr = node;
    } else {
        plotCurr->nxt = node;
        plotCurr = node;
    }
}

void pwlListAdd(circuitList *node, char *token) {
    node->pwl->t1 = (char *) malloc(sizeof (token));
    strcpy(node->pwl->t1, token);
    token = strtok(NULL, " \n\r()");
    //printf("PWL i: %s\n", token);
    node->pwl->i1 = (char *) malloc(sizeof (token));
    strcpy(node->pwl->i1, token);
    if (!listHasElement("pwl")) {
        pwlRoot = node->pwl;
        pwlCurr = node->pwl;
    } else {
        pwlCurr->nxt = node->pwl;
        pwlCurr = node->pwl;
    }
}

void populateList(char *list, char *token) {
    //printf("Populating %s\n", list);
    if (!strcmp(list, "circuit")) {
        int field = 0;
        int innerfield = 0;
        circuitList *l;
        l = (circuitList *) malloc(sizeof (circuitList));
        if (l != NULL) {
            l->transient_spec = "NO";
            while (token != NULL) {
                if (field == 0) {
                    l->name = (char *) malloc(sizeof (token));
                    strcpy(l->name, token);
                } else if (field == 1) {
                    l->positive = (char *) malloc(sizeof (token));
                    strcpy(l->positive, token);
                } else if (field == 2) {
                    l->negative = (char *) malloc(sizeof (token));
                    strcpy(l->negative, token);
                } else if (field == 3) {
                    l->value = (char *) malloc(sizeof (token));
                    strcpy(l->value, token);
                } else if (field == 4) {
                    //printf("list\n");
                    if (strstr(token, "EXP")) {
                        //printf("token: %s\n", token);
                        l->exp = (expList *) malloc(sizeof (expList));
                        l->transient_spec = "EXP";
                        while (token != NULL) {
                            if (innerfield == 1) {
                                //printf("innerfield: %d\n", innerfield);
                                l->exp->i1 = (char *) malloc(sizeof (token));
                                //printf("token: %s\n", token);
                                strcpy(l->exp->i1, token);
                            } else if (innerfield == 2) {
                                //printf("innerfield: %d\n", innerfield);
                                l->exp->i2 = (char *) malloc(sizeof (token));
                                //printf("token: %s\n", token);
                                strcpy(l->exp->i2, token);
                            } else if (innerfield == 3) {
                                //printf("innerfield: %d\n", innerfield);
                                l->exp->td1 = (char *) malloc(sizeof (token));
                                //printf("token: %s\n", token);
                                strcpy(l->exp->td1, token);
                            } else if (innerfield == 4) {
                                //printf("innerfield: %d\n", innerfield);
                                l->exp->tc1 = (char *) malloc(sizeof (token));
                                //printf("token: %s\n", token);
                                strcpy(l->exp->tc1, token);

                            } else if (innerfield == 5) {
                                //printf("innerfield: %d\n", innerfield);
                                l->exp->td2 = (char *) malloc(sizeof (token));
                                //printf("token: %s\n", token);
                                strcpy(l->exp->td2, token);

                            } else if (innerfield == 6) {
                                //printf("innerfield: %d\n", innerfield);
                                l->exp->tc2 = (char *) malloc(sizeof (token));
                                //printf("token: %s\n", token);
                                strcpy(l->exp->tc2, token);
                            }
                            //printf("exit\n");
                            innerfield++;
                            token = strtok(NULL, " \n\r()");
                            //printf("token: %s\n", token);
                        }
                    } else if (strstr(token, "SIN")) {
                        l->sin = (sinList *) malloc(sizeof (sinList));
                        //printf("token: %s\n", token);
                        l->transient_spec = "SIN";
                        innerfield = 0;
                        while (token != NULL) {
                            if (innerfield == 1) {
                                //printf("innerfield: %d\n", innerfield);
                                l->sin->i1 = (char *) malloc(sizeof (token));
                                strcpy(l->sin->i1, token);
                            } else if (innerfield == 2) {
                                //printf("innerfield: %d\n", innerfield);
                                l->sin->ia = (char *) malloc(sizeof (token));
                                strcpy(l->sin->ia, token);
                            } else if (innerfield == 3) {
                                //printf("innerfield: %d\n", innerfield);
                                l->sin->fr = (char *) malloc(sizeof (token));
                                strcpy(l->sin->fr, token);
                            } else if (innerfield == 4) {
                                //printf("innerfield: %d\n", innerfield);
                                l->sin->td = (char *) malloc(sizeof (token));
                                strcpy(l->sin->td, token);

                            } else if (innerfield == 5) {
                                //printf("innerfield: %d\n", innerfield);
                                l->sin->df = (char *) malloc(sizeof (token));
                                strcpy(l->sin->df, token);

                            } else if (innerfield == 6) {
                                //printf("innerfield: %d\n", innerfield);
                                l->sin->ph = (char *) malloc(sizeof (token));
                                strcpy(l->sin->ph, token);
                            }
                            innerfield++;
                            token = strtok(NULL, " \n\r()");
                        }
                    } else if (strstr(token, "PULSE")) {
                        l->pulse = (pulseList *) malloc(sizeof (pulseList));
                        //printf("token: %s\n", token);
                        l->transient_spec = "PULSE";
                        innerfield = 0;
                        while (token != NULL) {
                            if (innerfield == 1) {
                                //printf("innerfield: %d\n", innerfield);
                                l->pulse->i1 = (char *) malloc(sizeof (token));
                                strcpy(l->pulse->i1, token);
                            } else if (innerfield == 2) {
                                //printf("innerfield: %d\n", innerfield);
                                l->pulse->i2 = (char *) malloc(sizeof (token));
                                strcpy(l->pulse->i2, token);
                            } else if (innerfield == 3) {
                                //printf("innerfield: %d\n", innerfield);
                                l->pulse->td = (char *) malloc(sizeof (token));
                                strcpy(l->pulse->td, token);
                            } else if (innerfield == 4) {
                                //printf("innerfield: %d\n", innerfield);
                                l->pulse->tr = (char *) malloc(sizeof (token));
                                strcpy(l->pulse->tr, token);
                            } else if (innerfield == 5) {
                                //printf("innerfield: %d\n", innerfield);
                                l->pulse->tf = (char *) malloc(sizeof (token));
                                strcpy(l->pulse->tf, token);
                            } else if (innerfield == 6) {
                                //printf("innerfield: %d\n", innerfield);
                                l->pulse->pw = (char *) malloc(sizeof (token));
                                strcpy(l->pulse->pw, token);
                            } else if (innerfield == 7) {
                                //printf("innerfield: %d\n", innerfield);
                                l->pulse->per = (char *) malloc(sizeof (token));
                                strcpy(l->pulse->per, token);
                            }
                            innerfield++;
                            token = strtok(NULL, " \n\r()");
                        }
                    } else if (strstr(token, "PWL")) {
                        l->pwl = (pwlList *) malloc(sizeof (pwlList));
                        l->transient_spec = "PWL";
                        token = strtok(NULL, " \n\r()");
                        while (token != NULL) {
                            //printf("PWL t: %s\n", token);
                            pwlListAdd(l, token);
                            token = strtok(NULL, " \n\r()");
                        }
                    }
                }
                field++;
                token = strtok(NULL, " \n\r()");
            }
            circuitListAdd(l);
        } else printf("Error: Node creation failed.\n");
    } else if (!strcmp(list, "mos")) {
        int field = 0;
        mosList *l;
        l = (mosList *) malloc(sizeof (mosList));
        if (l != NULL) {
            while (token != NULL) {
                if (field == 0) {
                    l->name = (char *) malloc(sizeof (token));
                    strcpy(l->name, token);
                } else if (field == 1) {
                    l->drain = (char *) malloc(sizeof (token));
                    strcpy(l->drain, token);
                } else if (field == 2) {
                    l->gain = (char *) malloc(sizeof (token));
                    strcpy(l->gain, token);
                } else if (field == 3) {
                    l->source = (char *) malloc(sizeof (token));
                    strcpy(l->source, token);
                } else if (field == 4) {
                    l->base = (char *) malloc(sizeof (token));
                    strcpy(l->base, token);
                } else if (field == 5) {
                    l->model = (char *) malloc(sizeof (token));
                    strcpy(l->model, token);
                }
                field++;
                token = strtok(NULL, " \n\r");
            }
            mosListAdd(l);
        } else printf("Error: Node creation failed.\n");
    } else if (!strcmp(list, "bjt")) {
        int field = 0;
        bjtList *l;
        l = (bjtList *) malloc(sizeof (bjtList));
        if (l != NULL) {
            while (token != NULL) {
                if (field == 0) {
                    l->name = (char *) malloc(sizeof (token));
                    strcpy(l->name, token);
                } else if (field == 1) {
                    l->collector = (char *) malloc(sizeof (token));
                    strcpy(l->collector, token);
                } else if (field == 2) {
                    l->base = (char *) malloc(sizeof (token));
                    strcpy(l->base, token);
                } else if (field == 3) {
                    l->emitter = (char *) malloc(sizeof (token));
                    strcpy(l->emitter, token);
                } else if (field == 4) {
                    l->model = (char *) malloc(sizeof (token));
                    strcpy(l->model, token);
                }
                field++;
                token = strtok(NULL, " \n\r");
            }
            bjtListAdd(l);
        } else printf("Error: Node creation failed.\n");
    } else if (!strcmp(list, "dc")) {
        hasDC = true;
        int field = 0;
        dc = (dcList *) malloc(sizeof (dcList));
        if (dc != NULL) {
            while (token != NULL) {
                if (field == 1) {
                    dc->value = (char *) malloc(sizeof (token));
                    strcpy(dc->value, token);
                } else if (field == 2) {
                    dc->start = (char *) malloc(sizeof (token));
                    strcpy(dc->start, token);
                } else if (field == 3) {
                    dc->end = (char *) malloc(sizeof (token));
                    strcpy(dc->end, token);
                } else if (field == 4) {
                    dc->step = (char *) malloc(sizeof (token));
                    strcpy(dc->step, token);
                }
                field++;
                token = strtok(NULL, " \n\r");
            }
        } else printf("Error: Node creation failed.\n");
    } else if (!strcmp(list, "plot")) {
        todoPlot = true;
        int field = 0;
        while (token != NULL) {
            plotList *l = (plotList *) malloc(sizeof (plotList));
            if (l != NULL) {
                if (field > 0) {
                    l->value = (char *) malloc(sizeof (token));
                    strcpy(l->value, token);
                    plotListAdd(l);
                }
                field++;
                token = strtok(NULL, " \n\riIvV()");
            } else printf("Error: Node creation failed.\n");
        }
    } else if (!strcmp(list, "tran")) {
        //printf("tran\n");
        todoTransient = true;
        int field = 0;
        tran = (tranList *) malloc(sizeof (tranList));
        if (tran != NULL) {
            while (token != NULL) {
                if (field == 1) {
                    tran->step = (char *) malloc(sizeof (token));
                    strcpy(tran->step, token);
                } else if (field == 2) {
                    tran->end = (char *) malloc(sizeof (token));
                    //printf("token: %s\n", token);
                    strcpy(tran->end, token);
                }
                field++;
                token = strtok(NULL, " \n\r");
            }
        } else printf("Error: Node creation failed.\n");
    }
}
