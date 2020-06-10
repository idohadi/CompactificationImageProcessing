// TODO: write docs

/** 
 * Calculate the Clebsch-Gordan coefficients of order (l1, l2, l, m)
 * 
 * MATLAB call form:
 *      output = cgVector(l1, l2, l, m)
 *  where
 *      rotation is a quaternion representaiton of a rotation R and vectors is an 3 x N MATLAB array.
 *      output is a 3 x N MATLAB array such that
 *          output(:, n) = R * vectors(:, n)
 * 
 * NOTE:
 *  Code performs no input checks.
  */

#include <stdint.h>
#include "mex.h"
#include "clebsch_gordan_coefficients.h"


long l1, l2, l, m;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Obtain input
    l1 = mxGetScalar(prhs[0]);
    l2 = mxGetScalar(prhs[1]);
    l = mxGetScalar(prhs[2]);
    m = mxGetScalar(prhs[3]);

    // Generate output
    plhs[0] = mxCreateDoubleMatrix(clebsch_gordan_upper_bound(l1, l2, m) - clebsch_gordan_lower_bound(l1, l2, m) + 1, 1, mxREAL);
    cg_vector(l1, l2, l, m, mxGetDoubles(plhs[0]));
}
