char *toString(int integer);
float calcItol(char *line);
char *identifyElement(char *name);
int identifyCircuit(char *name);
void printVector(char *name, double *vector, size_t size);
void printPlot(plotList *l, double *vector);
float gslNorm(gsl_vector *vector);
float norm(double *vector, size_t size);
double dotProduct(double *a, double *b, size_t size);
double *addVectors(double *a, double *b, size_t size);
double *mulConstant(double *a, double c, size_t size);
cs *composeDiagonal(cs *A);
void printGslMatrix(char *name, gsl_matrix *matrix);
void printGslVector(char *name, gsl_vector *vector);
void printGslPlot(plotList *l, gsl_vector *vector);
double *solvePreconditioner(cs *M, double *z, double * r);
gsl_matrix *composeGslDiagonal(gsl_matrix *A);
gsl_vector *solveGslPreconditioner(gsl_matrix *M, gsl_vector *z, gsl_vector * r);
gsl_vector *computeGslAp(gsl_matrix *A, gsl_vector *p);
double *gslVectorToDouble(gsl_vector *gVector);
gsl_vector *doubleToGslVector(double *dVector, size_t size);
float pulse_solve(int k, float t, pulseList * node);
float exp_solve(float t, expList * node);
float sin_solve(float t, sinList * node);
float pwl_solve(float t, pwlList * node);
void LU_solve(gsl_matrix *matrixA, gsl_vector *vectorX, gsl_vector *vectorB);