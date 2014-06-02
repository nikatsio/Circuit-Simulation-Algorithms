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

#define CHUNK 1024
#define _XOPEN_SOURCE 500

bool todoCholesky = false;
bool todoCG = false;
bool todoBiCG = false;
bool todoCSCholesky = false;
bool todoCSCG = false;
bool todoCSBiCG = false;
bool todoCSLU = false;
bool todoTransient = false;
bool todoBE = false;
bool todoPlot = false;
float itol;

void parseFile(char url[]) {
    char *line = malloc(CHUNK);
    if (line != NULL) {
        FILE *file = fopen(url, "r");
        if (file) {
            while (fgets(line, CHUNK, file)) {
                if (line[0] != '*' && line[0] != '\0' && line[0] != '\r' && line[0] != '\n' && line[0] != '\t' && line[0] != ' ') {
                    if (strstr(line, ".OPTIONS SPARSE ITER") != NULL || strstr(line, ".options sparse iter") != NULL) {
                        if (strstr(line, ".OPTIONS SPARSE ITER SPD") != NULL || strstr(line, ".options sparse iter spd") != NULL) todoCSCG = true;
                        else todoCSBiCG = true;
                        itol = calcItol(line);
                    } else if (strstr(line, ".OPTIONS SPARSE SPD") != NULL || strstr(line, ".options sparse spd") != NULL) todoCSCholesky = true;
                    else if (strstr(line, ".OPTIONS SPARSE") != NULL || strstr(line, ".options sparse") != NULL) todoCSLU = true;
                    else if (strstr(line, ".OPTIONS ITER") != NULL || strstr(line, ".options iter") != NULL) {
                        if (strstr(line, ".OPTIONS ITER SPD") != NULL || strstr(line, ".options iter spd") != NULL) todoCG = true;
                        else todoBiCG = true;
                        itol = calcItol(line);
                    } 
                    else if (strstr(line, ".OPTIONS SPD") != NULL || strstr(line, ".options spd") != NULL) todoCholesky = true;
                    else if (strstr(line, ".OPTIONS METHOD=BE") != NULL || strstr(line, ".options method=be") != NULL) todoBE = true;
                    else if (strstr(line, ".OPTIONS METHOD=TR") != NULL || strstr(line, ".options method=tr") != NULL) {}
                    else {
                        char *token = strtok(line, " \n\r");
                        populateList(identifyElement(token), token);
                    }
                }
            }
            if (ferror(file)) printf("Error: ferror().\n");
            fclose(file);
        } else printf("Error: File not found.\n");
    } else printf("Error: malloc() failure.\n");
}