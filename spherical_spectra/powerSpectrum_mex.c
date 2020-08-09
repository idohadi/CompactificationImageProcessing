// TODO: plan the new version

/** 
 * Calculate the bispectrum and its gradient.
 * 
 * MATLAB call form:
 *      p = powerSpectrum_mex(shc, bandlimit)
 *      p = powerSpectrum_mex(shc, bandlimit)
 *  where
 *      shc         column or row complex array of length (bandlimit+1)^2 of spherical harmonics coefficients
 *      bandlimit   scalar, the bandlimit of the function represented by shc
 *      p           bandlimit+1 column array, the power spectrum of shc
 * 
 * NOTES:
 *  (1) The code performs no input checks.
 *  (2) This function is wrapped by powerSpectrum.m.
 * 
  */

#include <stdint.h>
#include "mex.h"
#include "clebsch_gordan_coefficients.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Obtain input
    // TODO

    // Handle first run
    // TODO

    // Create output
    // TODO

    // Compute output
    // TODO
}
