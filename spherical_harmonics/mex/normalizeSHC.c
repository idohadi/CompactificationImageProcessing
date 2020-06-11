// TODO: docs in doxygen
/** 
 * Normalizes an array of spherical harmonics coefficient so that its power spctrum is all ones.
 * 
 * MATLAB call form:
 *      nshc = normalizeSHC(shc, bandlimit, real_vs_imag)
 *  where
 *      shc is a column array of size (bandlimit+1)^2 if real_vs_imag=0 and of size 2*(bandlimit+1)^2 if real_vs_imag~=0.
 *      powspec is always a column array of size bandlimit+1.
 * 
 * NOTE:
 *  Code performs no input checks.
  */

#include <stdint.h>
#include <stdbool.h>
#include "mex.h"
#include "spherical_harmonics.h"

size_t bandlimit;
size_t real_vs_imag;

r_shc shc1;
r_shc output_shc1;
c_shc shc2;
c_shc output_shc2;


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Get input variables
    bandlimit = (size_t) mxGetScalar(prhs[1]);
    real_vs_imag = (size_t) mxGetScalar(prhs[2]);

    // Create output
    if (real_vs_imag==0)
    {
        plhs[0] = mxCreateDoubleMatrix((bandlimit+1)*(bandlimit+1), 1, mxREAL);
        shc1 = r_init_shc(bandlimit, mxGetDoubles(prhs[0]));
        output_shc1 = r_init_shc(bandlimit, mxGetDoubles(plhs[0]));
        r_normalize_shc(&shc1, &output_shc1);
    }
    else
    {
        plhs[0] = mxCreateDoubleMatrix(2*(bandlimit+1)*(bandlimit+1), 1, mxREAL);
        shc2 = c_init_shc(bandlimit, mxGetDoubles(prhs[0]));
        output_shc2 = c_init_shc(bandlimit, mxGetDoubles(plhs[0]));
        c_normalize_shc(&shc2, &output_shc2);
    }
}
