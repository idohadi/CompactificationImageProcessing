// TODO: write docs in doxygen

#include "alegendre.h"

void alegendre_init()
{
    double dsize;
    alegendre_eval_init_wrapper(&dsize);
}

double alegendre(const long l, const long m, const double x)
{
    double alpha, alphader, vallogp, vallogq, valp, valq;
    alegendre_eval_wrapper(l, -m, x, 
                            &alpha, &alphader, 
                            &vallogp, &vallogq, 
                            &valp, &valq);
    return valp;
}

long alegendre_root_count(const long l, const long m)
{
    long nproots, nqroots;
    alegendre_nroots_wrapper(l, -m, &nproots, &nqroots);
    return nproots;
}

double alegendre_root(const long l, const long m, const long root_order)
{
    double x;
    alegendre_proot_wrapper(l, -m, root_order, &x);
    return x;
}
