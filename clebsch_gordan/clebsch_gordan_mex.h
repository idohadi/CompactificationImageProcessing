#include "mex.h"


typedef double ***** CGTable;


void create_CGTable(CGTable * const cgt, const mxArray ** CGs, const size_t bl);

void destroy_CGTable(CGTable * const cgt, const size_t bl);

double get_cg(CGTable * const cgt, const long l1, const long l2, const long l, const long m, const long m1);
