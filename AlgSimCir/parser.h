extern bool todoCholesky;
extern bool todoCG;
extern bool todoBiCG;
extern bool todoCSCholesky;
extern bool todoCSCG;
extern bool todoCSBiCG;
extern bool todoCSLU;
extern bool todoTransient;
extern bool todoPlot;
extern bool todoBE;
extern float itol;

void parseLine(char line[]);
void parseFile(char url[]);