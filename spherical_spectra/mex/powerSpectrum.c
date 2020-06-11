// TODO: docs in doxygen
/** 
 * Return the power spectrum of a bandlimited spherical function represnted by its spherical harmonics coefficients. 
 * 
 * MATLAB call form:
 *      powspec = powerSpectrum(shc, bandlimit, real_vs_imag)
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
#include "spherical_spectra.h"

size_t bandlimit;
size_t real_vs_imag;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Get input variables
    bandlimit = (size_t) mxGetScalar(prhs[1]);
    real_vs_imag = (size_t) mxGetScalar(prhs[2]);

    // Create output
    plhs[0] = mxCreateDoubleMatrix(bandlimit+1, 1, mxREAL);
    if (real_vs_imag==0)
    {
        r_shc shc = r_init_shc(bandlimit, mxGetDoubles(prhs[0]));
        r_power_spectrum(&shc, mxGetDoubles(plhs[0]));
    }
    else
    {
        c_shc shc = c_init_shc(bandlimit, mxGetDoubles(prhs[0]));
        c_power_spectrum(&shc, mxGetDoubles(plhs[0]));
    }
}