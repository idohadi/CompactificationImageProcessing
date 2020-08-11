/** 
 * Normalize spherical harmonics coefficients array to have a power spectrum of all ones.
 * 
 * MATLAB call form:
 *      nshc = normSHC_mex(shc, shcAbsSqr, bandlimit)
 *  where
 *      shc         (bandlimit+1)^2 x N complex array such that shc(:, n) are spherical 
 *                  harmoncis coefficients of a function of bandlimit
 *      shcAbsSqr   (bandlimit+1)^2 x N complex array such that 
 *                      shcAbsSqr(j, n) = shc(j, n) * conj(shc(j, n))
 *      bandlimit   scalar, the bandlimit of the function represented by the spherical harmonics cofficients
 * 
 * NOTES:
 *  (1) The code performs no input checks.
 *  (2) This function is wrapped by normSHC.m.
 * 
  */

#include <stdint.h>
#include <math.h>
#include "mex.h"


mxComplexDouble *shc;
double *shcAbsSqr;
size_t bandlimit;

mxComplexDouble *nshc;

long N;
long M;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Obtain input
    shc = mxGetComplexDoubles(prhs[0]);
    shcAbsSqr = mxGetDoubles(prhs[1]);
    bandlimit = mxGetScalar(prhs[2]);
    
    M = mxGetM(prhs[0]);
    N = mxGetN(prhs[0]);

    // Create output
    plhs[0] = mxCreateDoubleMatrix(M, N, mxCOMPLEX);
    nshc = mxGetComplexDoubles(plhs[0]);
    for (long n = 0; n<N*M; ++n)
    {
        nshc[n] = shc[n];
    }
    
    // Compute output
    double norm;
    for (size_t n = 0; n<N; ++n)
    {
        for (long l = 0; l<=bandlimit; ++l)
        {
            norm = 0;
            for (long m = -l; m<=l; ++m)
            {
                norm += shcAbsSqr[n*M + l*(l+1) + m];
            }
            norm = sqrt(norm);

            for (long m = -l; m<=l; ++m)
            {
                nshc[n*M + l*(l+1) + m].real /= norm;
                nshc[n*M + l*(l+1) + m].imag /= norm;
            }
        }
    }
}
 