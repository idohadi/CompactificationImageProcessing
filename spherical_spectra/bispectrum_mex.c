// TODO: plan the new version

/** 
 * Calculate the bispectrum and its gradient.
 * 
 * MATLAB call form:
 *      b = bispectrum_mex(shc, bandlimit, CGs)
 *      [b, grad] = bispectrum_mex(shc, bandlimit, CGs)
 *  where
 *      shc         column or row complex array of length (bandlimit+1)^2 of spherical harmonics coefficients
 *      bandlimit   scalar, the bandlimit of the function represented by shc
 *      CGs         a cell array containing the precomputed Clebsch-Gordan coefficents for bandlimit
 *      b           bispectrum vector, containing the bispectrum invariants b_{l1,l2,l} for 
 *                      0<=l1<=bandlimit, 0<=l2<=l1, abs(l1-l2)<=l<=min(bandlmit, l1+l2)
 *                  ordered in lexigraphical order.
 *      grad        the transpose of the gradient of the mapping shc->b
 * 
 * NOTES:
 *  (1) The code performs no input checks.
 *  (2) This function is wrapped by bispectrum.m.
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
