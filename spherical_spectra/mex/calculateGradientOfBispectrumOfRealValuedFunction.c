// TODO
/** 
 * Return the gradient of the bispectrum vector based on spherical harmonics vector. 
 * 
 * MATLAB call form:
 *      bispectrum_gradient = calculateGradientOfBispectrumOfRealValuedFunction(shc, bandlimit)
 *  where
 *      bispectrum_gradient is the gradient (its transpose, actually).
 *      
 *      shc is a column or row vector of spherical harmonics coefficients of a real-valued spherical function
 *          shc 
 *              =   [f_{r, 0, 0}, 
 *                  f_{r, 1, 0}, f_{r, 1, 1}, f_{i, 1, 1}, 
 *                  f_{r, 2, 0}, f_{r, 2, 1}, f_{i, 2, 1}, f_{r, 2, 2}, f_{i, 2, 2}, ...
 *                  f_{r, bandlimit, 0}, f_{r, bandlimit, 1}, f_{i, bandlimit, 1}, ..., f_{r, bandlimit, bandlimit}, f_{i, bandlimit, bandlimit} ]
 *          f_{r,l,m} is the real part of the spherical harmonics coefficient of order l, degree m
 *          f_{i,l,m} is its imaginary part.
 * 
 * NOTE:
 *  Code performs no input checks.
  */


#include "mex.h"
#include "..\c\clebsch_gordan_coefficients.h"
#include "..\c\spherical_spectra.h"
#include <stdint.h>

bool first_run = true;
size_t ***bispectrum_lookup_table;
cb_table cb_coeffs;
size_t bandlimit = 0;
size_t previous_bandlimit = 0;
double *shc;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Get input variables
    shc = mxGetDoubles(prhs[0]);
    bandlimit = (size_t) mxGetScalar(prhs[1]);

    // Handle first run
    if (first_run == true)
    {
        mexPrintf("First run. Initializing Clebsch-Gordan coefficients and bispectrum lookup table.\n");
        previous_bandlimit = bandlimit;
        cb_coeffs = allocate_memeory_for_clebsch_gordan_coefficients(bandlimit);
        calculate_clebsch_gordan_coefficients(bandlimit, cb_coeffs);
        bispectrum_lookup_table = r_build_bispectrum_lookup_table(bandlimit);
        mexPrintf("Initialization completed.\n");
        first_run = false;
    }
    // Handle chagne in bandlimit between runs
    if (previous_bandlimit!=bandlimit)
    {
        mexPrintf("Bandlimit chagned. Reinitializing Clebsch-Gordan coefficients and bispectrum lookup table.\n");
        free_memory_for_clebsch_gordan_coefficients(previous_bandlimit, cb_coeffs);
        r_destroy_bispectrum_lookup_table(bispectrum_lookup_table, previous_bandlimit);

        previous_bandlimit = bandlimit;
        cb_coeffs = allocate_memeory_for_clebsch_gordan_coefficients(bandlimit);
        calculate_clebsch_gordan_coefficients(bandlimit, cb_coeffs);
        bispectrum_lookup_table = r_build_bispectrum_lookup_table(bandlimit);
        mexPrintf("Reinitialization completed.\n");
    }

    // Create output
    plhs[0] = mxCreateDoubleMatrix((bandlimit+1)*(bandlimit+1), bispectrum_lookup_table[bandlimit][bandlimit][bandlimit]+1, mxREAL);

    // Calculate the bispectrm / its gradient
    r_calculate_bispectrum_gradient(shc, bandlimit, bispectrum_lookup_table, cb_coeffs, mxGetDoubles(plhs[0]));
}