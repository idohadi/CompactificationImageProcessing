#include <stdlib.h>

/* Typedef */
typedef double***** clebsch_gordan_table;
typedef 
    struct cg_table
    {
        clebsch_gordan_table table;
        size_t bandlimit;
    }  cg_table;



/* Function declarations */
long cg_lower_bound(const long l1, const long l2, const long m);


long cg_upper_bound(const long l1, const long l2, const long m);


double backward_substitution_central_diagonal(const long l1, const long l2, const long l, const long m1, const long m);


double backward_substitution_secondary_diagonal(const long l1, const long l2, const long m1, const long m);


void cg_vector(const long l1, const long l2, const long l, const long m, double * restrict cg);


double get_cg(const cg_table *table, const long l1, const long l2, const long l, const long m, const long m1);
