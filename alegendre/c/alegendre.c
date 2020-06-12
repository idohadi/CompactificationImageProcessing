// TODO: write docs in doxygen

#include <math.h>
#include <float.h>
#include "alegendre.h"

#define PI 3.14159265358979323846

void alegendre_init()
{
    double dsize;
    alegendre_eval_init_wrapper(&dsize);
}

double alegendre(const long l, const long m, const double t)
{
    /* 
    Calculates the associated Legendre polynomial P_{l}^{m} (cos(t)).

    Must run alegendre_init once before using this function.
    */

    double alpha, alphader, vallogp, vallogq, valp, valq;
    if (t>PI/2 || t<-PI/2)
    {
        return 0.0;
    }

    if (t>=0)
    {
        alegendre_eval_wrapper(l, -m, t, 
                            &alpha, &alphader, 
                            &vallogp, &vallogq, 
                            &valp, &valq);
    }
    else
    {
        alegendre_eval_wrapper(l, -m, -t, 
                            &alpha, &alphader, 
                            &vallogp, &vallogq, 
                            &valp, &valq);
    }
    
    if (valp<=DBL_EPSILON && valp>=-DBL_EPSILON)
    {
        return 0.0;
    }
    else
    {
        return valp/sqrt(2.0*PI*sin(t));
    }
    
    return valp;
}

long alegendre_root_count(const long l, const long m)
{
    /* 
    TODO: docs

    Must run alegendre_init once before using this function.
    */
    long nproots, nqroots;
    alegendre_nroots_wrapper(l, -m, &nproots, &nqroots);
    return nproots;
}

double alegendre_root(const long l, const long m, const long root_order)
{
    /* 
    TODO: document and also make sure to understand exactly what it does and that the return value is what I want.

    Must run alegendre_init once before using this function.
    */
    double x;
    alegendre_proot_wrapper(l, -m, root_order, &x);
    return x;
}
