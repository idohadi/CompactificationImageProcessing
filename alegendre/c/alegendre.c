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

double alegendre(const long l, const long m, double t)
{
    /* 
    Calculates the normalized associated Legendre polynomial 
        S_{l}^{m} (t) := sqrt(((2l+1)/(4*pi)) * ((l-m)! / (l+m)!))    *   P_{l}^{m} (cos(t)) 
    for 0 <= t <= pi.
    For t>pi/2, it returns the value at pi/2.
    For t<-pi/2, it returns the value at -pi/2.

    NOTES:
        (1) Must run alegendre_init once before using this function.
        (2) See Bateman, 1953, 3.4(17), 3.4(20), 3.9.2(9), 3.9.2(13)-(15)
    */

    double alpha, alphader, vallogp, vallogq, valp, valq;

    // Find t in [0,pi] such that cos(t_original)=cos(t_new)
    // Normalize t to be in [0,2pi]
    if (t>2.0*PI || t<0.0)
    {
        t = (t/(2.0*PI) - floor((t/(2.0*PI))))*2.0*PI;
    }
    // Normalize t to be in [0,pi]
    if (t>PI)
    {
        // t = 2*PI-t, in place
        t *= -1;
        t += 2.0*PI;
    }
    
    // Calculate S_{l}^{-|m|} (t)
    if (t>0.0 && t<=0.5*PI)
    {
        // Calcualte S_{l}^{-|m|}(t) * sqrt(sin(t))
        alegendre_eval_wrapper(l, (m>=0 ? m : -m), t, 
                            &alpha, &alphader, 
                            &vallogp, &vallogq, 
                            &valp, &valq);
        valp /= sqrt(2.0*PI*sin(t));
    }
    else if (t>0.5*PI && t<PI)
    {
        // Calcualte S_{l}^{-|m|}(t-pi/2) * sqrt(sin(t-pi/2))
        alegendre_eval_wrapper(l, (m>=0 ? m : -m), t-0.5*PI, 
                            &alpha, &alphader, 
                            &vallogp, &vallogq, 
                            &valp, &valq);
        valp /= sqrt(2.0*PI*sin(t-0.5*PI));
        // Calculate S_{l}^{-|m|} (t); see Bateman 3.4(19)
        valp *= (((l+m) % 2) == 0 ? 1 : -1);
    }
    else if (t==0)
    {
        // Calculate S_{l}^{-|m|}(0) and return it; see Bateman, 3.9.2(9)
        if (m==0)
        {
            return sqrt((2.0*l+1.0)/(4*PI));
        }
        else
        {
            return 0.0;
        }
    }
    else // if (t == PI), the only remaining case
    {
        // Calculate S_{l}^{-|m|}(PI) and return it; see Bateman, 3.9.2(9) and 3.4(19)
        if (m==0)
        {
            return ((m%2)==0 ? 1 : -1)*sqrt((2.0*l+1.0)/(4*PI));
        }
        else
        {
            return 0.0;
        }
    }
    
    if (m<=0)
    {
        return valp;
    }
    else
    {
        // Calcualte S_{l}^{m}(t) for m>0; see Bateman, 3.4(17)
        valq = gamma_ratio(l, m);
        valp *= ((m%2)==0 ? 1 : -1) * valq*valq;
        return valp;
    }
    
    return valp;
}


double alegendre2(const long l, const long m, const double x)
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