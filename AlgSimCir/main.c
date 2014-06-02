#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdbool.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include "csparse.h"
#include "parser.h"
#include "list.h"
#include "hashtbl.h"
#include "table.h"
#include "util.h"

int main(int argc, char *argv[]) {
    listInit("circuit");
    listInit("mos");
    listInit("bjt");
    listInit("terminal");
    parseFile(argv[1]);
    if (listHasElement("circuit")) {
        hashtbl = populateHashtbl(circuitRoot);
        if (todoTransient) {
            allocateTransientSystem();
            populateTransientSystem(circuitRoot, hashtbl, tilda_G, tilda_C, b);
            solveTransientSystem(circuitRoot);
        } else {
            if (todoCSCG || todoCSBiCG || todoCSCholesky || todoCSLU) {
                calcNonZeroElem(circuitRoot, hashtbl);
                allocateCSSystem();
                populateCSSystem(circuitRoot, hashtbl, cs_A);
                solveCSSystem();
            } else {
                allocateSystem();
                populateSystem(circuitRoot, hashtbl, A, b);
                solveSystem();
            }
        }/**/
        return 0;
    }
}