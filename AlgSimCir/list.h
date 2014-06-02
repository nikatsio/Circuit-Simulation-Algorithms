typedef struct circuitNode {
    char *name;
    char *positive;
    char *negative;
    char *value;
    char *transient_spec;
    struct expNode *exp;
    struct sinNode *sin;
    struct pulseNode *pulse;
    struct pwlNode *pwl;
    struct circuitNode *nxt;
} circuitList;

typedef struct mosNode {
    char *name;
    char *drain;
    char *gain;
    char *source;
    char *base;
    char *model;
    struct mosNode *nxt;
} mosList;

typedef struct bjtNode {
    char *name;
    char *collector;
    char *base;
    char *emitter;
    char *model;
    struct bjtNode *nxt;
} bjtList;

typedef struct dcNode {
    char *value;
    char *start;
    char *end;
    char *step;
} dcList;

typedef struct plotNode {
    char *value;
    struct plotNode *nxt;
} plotList;

typedef struct tranNode {
    char *step;
    char *end;
} tranList;

typedef struct expNode {
    char *i1;
    char *i2;
    char *td1;
    char *tc1;
    char *td2;
    char *tc2;
} expList;

typedef struct sinNode {
    char *i1;
    char *ia;
    char *fr;
    char *td;
    char *df;
    char *ph;
} sinList;

typedef struct pulseNode {
    char *i1;
    char *i2;
    char *td;
    char *tr;
    char *tf;
    char *pw;
    char *per;
} pulseList;

typedef struct pwlNode {
    char *t1;
    char *i1;
    struct pwlNode *nxt;
} pwlList;

extern circuitList *circuitRoot;
extern circuitList *circuitCurr;
extern mosList *mosRoot;
extern mosList *mosCurr;
extern bjtList *bjtRoot;
extern bjtList *bjtCurr;
extern dcList *dc;
extern tranList *tran;
extern plotList *plotRoot;
extern plotList *plotCurr;
extern pwlList *pwlRoot;
extern pwlList *pwlCurr;

extern int element;
extern bool hasDC;

void listInit(char *list);
int listHasElement(char *list);
void circuitListAdd(circuitList *node);
void *circuitListGet(char *name);
void mosListAdd(mosList *node);
void bjtListAdd(bjtList *node);
void plotListAdd(plotList *node);
void pwlListAdd(circuitList *node, char *token);
void populateList(char *list, char *token);