// TODO
/** 
 * Return an array of random normalized spherical harmonics coefficients of a given bandlimit. 
 * 
 * MATLAB call form:
 *      shc = randomNormalizedSHC(bandlimit, real_vs_imag)
 *  where
 *      shc is a column array of size (bandlimit+1)^2 if real_vs_imag=0 and of size 2*(bandlimit+1)^2 if real_vs_imag~=0
 * 
 * NOTE:
 *  Code performs no input checks.
  */

#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include "mex.h"
#include "spherical_harmonics.h"
#include "SFMT.h"


size_t bandlimit;
size_t real_vs_imag;
sfmt_t sfmt;
bool first_run = true;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Initilize SFMT on the first run
    if (first_run==true)
    {
        sfmt_init_gen_rand(&sfmt, time(NULL));
        first_run = false;
    }

    // Get input variables
    bandlimit = (size_t) mxGetScalar(prhs[0]);
    real_vs_imag = (size_t) mxGetScalar(prhs[1]);

    // Create output
    if (real_vs_imag==0)
    {
        plhs[0] = mxCreateDoubleMatrix((bandlimit+1)*(bandlimit+1), 1, mxREAL);
        r_shc shc = r_init_shc(bandlimit, mxGetDoubles(plhs[0]));
        r_random_normalized_shc(&sfmt, &shc);
    }
    else
    {
        plhs[0] = mxCreateDoubleMatrix(2*(bandlimit+1)*(bandlimit+1), 1, mxREAL);
        c_shc shc = c_init_shc(bandlimit, mxGetDoubles(plhs[0]));
        c_random_normalized_shc(&sfmt, &shc);
    }
}