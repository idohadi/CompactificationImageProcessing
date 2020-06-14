// TODO: write docs

/** 
 * Calculate the associated Legendre polynomial of input array of nubmers in [-1,1].
 * 
 * MATLAB call form:
 *      output = evalAlegendre(l, m, x)
 *  where
 *      l, m are order, degree of the polynomial.
 *      x is a double array. All numbers need to be in [-1, 1].
 *      output(j) is the associated Legendre polynomial if x(i) in [-1,1] and NaN otherwise.
 * 
 * NOTE:
 *  Code performs no input checks.
  */

#include <stdint.h>
#include <stdbool.h>
#include "mex.h"
#include "alegendre.h"


long l, m;
double *x; 
double *output;
size_t N; 
bool first_run = true;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // // Initialize the fast Legendre function
    if (first_run)
    {
        alegendre_init();
        first_run = false;
        mexPrintf("evalAlegendre initialized successfully.\n");
    }

    // Obtain input
    l = mxGetScalar(prhs[0]);
    m = mxGetScalar(prhs[1]);
    x = mxGetDoubles(prhs[2]);
    N = mxGetNumberOfElements(prhs[2]);

    // Generate output
    plhs[0] = mxCreateDoubleMatrix(N, 1, mxREAL);
    output = mxGetDoubles(plhs[0]);

    for (size_t i = 0; i<N; ++i)
    {
        if (x[i]<=1 && x[i]>=-1)
        {
            output[i] = alegendre2(l, m, x[i]);
        }
        else
        {
            output[i] = mxGetNaN();
        }
    }
}
