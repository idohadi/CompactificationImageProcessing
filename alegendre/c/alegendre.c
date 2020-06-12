// TODO: write docs in doxygen

#include <math.h>
#include <float.h>
#include <stdio.h>
#include "alegendre.h"

#define PI 3.14159265358979323846

void alegendre_init()
{
    double dsize;
    alegendre_eval_init_wrapper(&dsize);
}

double alegendre(const long l, const long m, double x)
{
    /* 
    Calculates the associated Legendre polynomial P_{l}^{m} (x).

    Must run alegendre_init once before using this function.
    */

    double alpha, alphader, vallogp, vallogq, valp, valq;

    if (!(x<1 && x>-1))
    {
        printf("You can evaluate the polynomial only on (-1,1).\n");
        return -100.0;
    }

    if (x<-1 || x>1)
    {
        return 0.0;
    }

    if (x>0 && x<1)
    {
        if (m<=0)
        {
            alegendre_eval_wrapper(l, -m, acos(x), 
                            &alpha, &alphader, 
                            &vallogp, &vallogq, 
                            &valp, &valq);
        }
        else
        {
            alegendre_eval_wrapper(l, m, acos(x), 
                            &alpha, &alphader, 
                            &vallogp, &vallogq, 
                            &valp, &valq);
            
            // TODO: doc why this is the formulas:
            valp *= ((m%2 == 0) ? 1 : -1)*gamma_ratio(l, m);
        }

        if (valp<=DBL_EPSILON && valp>=-DBL_EPSILON)
        {
            return 0.0;
        }
        else
        {
            return valp/sqrt(2.0*PI*sin(acos(x)));
        }
    }
    else if  (x>-1 && x<0)
    {
        if (m<=0)
        {
            alegendre_eval_wrapper(l, -m, acos(-x), 
                            &alpha, &alphader, 
                            &vallogp, &vallogq, 
                            &valp, &valq);
        }
        else
        {
            alegendre_eval_wrapper(l, m, acos(-x), 
                            &alpha, &alphader, 
                            &vallogp, &vallogq, 
                            &valp, &valq);
            
            // TODO: doc why this is the formulas:
            valp *= ((m%2 == 0) ? 1 : -1)*gamma_ratio(l, m);
        }

        valp *= ((l+m)%2 == 0 ? 1 : -1);
        if (valp<=DBL_EPSILON && valp>=-DBL_EPSILON)
        {
            return 0.0;
        }
        else
        {
            return valp/sqrt(2.0*PI*sin(acos(-x)));
        }
    }
    else
    {
        if ((l+m)%2 == 1)
        {
            return 0.0;
        }
        else
        {
            printf("For t = 0 and even l+m, this function is inaccuarte.\n");
            return -100.0;
        }
        
    }
    
    
    return -100.0;
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


double gamma_ratio(const long l, const long m)
{
    /* 
    TODO: docs

    Calculate Gamma(l+m+1) / Gamma(l-m+1)  = (l+m)!/(l-m)!
    */

   double val;
   alegendre_gamma_ratio2_wrapper(l, m, &val);
   return exp(val);
}

double double_factorial(const long l, const long m)
{
    /* 
    TODO: docs. It is for the calculation of the alegendre at 0

    Calculate (n+m-1)!! / (l-m+1)!!
    */

}