/** 
 * TODO: 
 *      (1) Update it to use the new code base for Clbesch-Gordan coeffs
 *      (2) Write better docs
 * 
 * MATLAB call form:
 *      V = buildV_mex(UUH, bandlimit)
 *  UUH = U*(U^*)
 * 
 * NOTE:
 *  Code performs no input checks.
  */


#include <stdint.h>
#include "mex.h"

size_t bandlimit = 0;
mxComplexDouble *UUH;

mxDouble *V;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Get input variables
    UUH = mxGetComplexDoubles(prhs[0]);
    bandlimit = (size_t) mxGetScalar(prhs[1]);
    long sz = mxGetN(prhs[0]);

    // Generate output
    plhs[0] = mxCreateDoubleMatrix(bandlimit+1, 1, mxREAL);
    V = mxGetDoubles(plhs[0]);

    // Build V
    for (long l = 0; l<=bandlimit; ++l)
    {
        for (long m = -l; m<=l; ++m)
        {
            V[l] += UUH[sz*(l*(l+1) + m) + l*(l+1) + m].real;
        }
        
    }
}
