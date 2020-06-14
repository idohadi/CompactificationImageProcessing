// TODO

/* Fortran functions. */
void alegendre_eval_init_wrapper(double *dsize);

void alegendre_eval_wrapper(double dnu, double dmu, double t, 
                            double * restrict alpha, double * restrict alphader, 
                            double * restrict vallogp, double * restrict vallogq, 
                            double * restrict valp, double * restrict valq);

void alegendre_nroots_wrapper(double dnu, double dmu, long *nproots, long *nqroots);

void alegendre_proot_wrapper(double dnu, double dmu, long j, double *t);

void alegendre_gamma_ratio2_wrapper(double dnu, double dmu, double * restrict val);

/* My C functions */
void alegendre_init();

double alegendre(const long l, const long m, double t);

double alegendre2(const long l, const long m, const double x);

long alegendre_root_count(const long l, const long m);

double alegendre_root(const long l, const long m, const long root_order);

double gamma_ratio(const long l, const long m);

double double_factorial(const long l, const long m);
