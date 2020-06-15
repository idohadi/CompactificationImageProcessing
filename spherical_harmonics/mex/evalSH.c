// TODO: write docs in doxygen
/** 
 * Evaluate real-valued spherical function, represented by spherical harmonics coefficients.
 * 
 * MATLAB call form:
 *      values = evalSHC(shc, bandlimit, theta, phi, real_vs_imag)
 *  where
 *      shc is a column array of size (bandlimit+1)^2 if real_vs_imag=0 and of size 2*(bandlimit+1)^2 if real_vs_imag~=0
 *      theta and phi and arrays of size N x 1 or 1 x N
 *      real_vs_imag is 0, if shc represents real-valued spherical function and non-zero if it represents a complex-valued spherical function.
 * 
 * NOTE:
 *  Code performs no input checks.
  */

#include <stdint.h>
#include <stdbool.h>
#include "mex.h"
#include "spherical_harmonics.h"
#include "alegendre.h"

r_shc shc;
double *theta;
double *phi;

size_t bandlimit;
size_t real_vs_imag;

size_t N;

double *values; 

bool first_run = true;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (first_run)
    {
        alegendre_init();
        first_run = false;
        mexPrintf("Algendere succesfully initialized.\n");
    }

    // Get input variables
    bandlimit = (size_t) mxGetScalar(prhs[1]);
    theta = mxGetDoubles(prhs[2]);
    phi = mxGetDoubles(prhs[3]);
    real_vs_imag = (size_t) mxGetScalar(prhs[4]);

    N = mxGetNumberOfElements(prhs[2])>=mxGetNumberOfElements(prhs[3]) 
            ? mxGetNumberOfElements(prhs[3]) 
            : mxGetNumberOfElements(prhs[2]);

    // Create output
    if (real_vs_imag==0)
    {
        plhs[0] = mxCreateDoubleMatrix(N, 1, mxREAL);
        values = mxGetDoubles(plhs[0]);

        shc = r_init_shc(bandlimit, mxGetDoubles(prhs[0]));
        for (size_t i = 0; i<N; ++i)
        {
            r_eval_sf(&shc, theta[i], phi[i], &values[i]);
        }
    }
    else
    {
        // TODO: handle complex-valued functions in future
    }
}