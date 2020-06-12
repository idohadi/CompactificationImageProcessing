// TODO

/* Fortran functions. */
void alegendre_eval_init_wrapper(double *dsize);

void alegendre_eval_wrapper(double dnu, double dmu, double t, 
                            double *alpha, double *alphader, 
                            double *vallogp, double *vallogq, 
                            double *valp, double *valq);

void alegendre_nroots_wrapper(double dnu, double dmu, long *nproots, long *nqroots);

void alegendre_proot_wrapper(double dnu, double dmu, long j, double *t);


/* My C functions */
void alegendre_init();

double alegendre(const long l, const long m, const double x);

long alegendre_root_count(const long l, const long m);

double alegendre_root(const long l, const long m, const long root_order);
