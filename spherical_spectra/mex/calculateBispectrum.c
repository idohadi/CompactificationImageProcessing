// TODO: docs in doxygen
/** 
 * Return a bispectrum vector based on spherical harmonics vector. 
 * 
 * MATLAB call form:
 *      [bispectrum, bispectrum_gradient] = calculateBispectrum(shc, bandlimit, real_vs_imag)
 *  where
 *      bispectrum is the bispectrum vector. 
 *      bispectrum_gradient (optional argument) is its gradient (its transpose, actually).
 *      shc is assumed to represent real-valued spherical function if real_vs_imag=0 and complex-valued otherwise.
 * 
 * NOTE:
 *  Code performs no input checks.
  */


#include <stdint.h>
#include "mex.h"
#include "spherical_spectra.h"

bool first_run = true;
r_blt r_lookup;
c_blt c_lookup;
cg_table cgs;
size_t bandlimit = 0;
size_t previous_bandlimit = 0;
size_t real_vs_imag;

r_shc shc1;
c_shc shc2;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Get input variables
    bandlimit = (size_t) mxGetScalar(prhs[1]);

    // Handle first run
    if (first_run == true)
    {
        mexPrintf("First run. Initializing Clebsch-Gordan coefficients and bispectrum lookup table.\n");
        previous_bandlimit = bandlimit;
        cgs = allocate_cg_table(bandlimit);
        calculate_cg_table(bandlimit, &cgs);
        r_lookup = r_build_bispectrum_lookup_table(bandlimit);
        c_lookup = c_build_bispectrum_lookup_table(bandlimit);
        mexPrintf("Initialization completed.\n");
        first_run = false;
    }

    // Handle chagne in bandlimit between runs
    if (previous_bandlimit!=bandlimit)
    {
        mexPrintf("Bandlimit chagned. Reinitializing Clebsch-Gordan coefficients and bispectrum lookup table.\n");
        // Destroying the lookup tables and Clebsch-Gordan table
        destroy_cg_table(&cgs);
        r_destroy_bispectrum_lookup_table(r_lookup, previous_bandlimit);
        c_destroy_bispectrum_lookup_table(c_lookup, previous_bandlimit);

        // Generate new tables
        previous_bandlimit = bandlimit;
        cgs = allocate_cg_table(bandlimit);

        calculate_cg_table(bandlimit, &cgs);
        r_lookup = r_build_bispectrum_lookup_table(bandlimit);
        c_lookup = c_build_bispectrum_lookup_table(bandlimit);

        mexPrintf("Reinitialization completed.\n");
    }

    // Create output
    if (real_vs_imag==0)
    {
        plhs[0] = mxCreateDoubleMatrix(r_lookup[bandlimit][bandlimit][bandlimit]+1, 1, mxREAL);
        shc1 = r_init_shc(bandlimit, mxGetDoubles(prhs[0]));

        if (nlhs>1)
        {
            plhs[1] = mxCreateDoubleMatrix((bandlimit+1)*(bandlimit+1), r_lookup[bandlimit][bandlimit][bandlimit]+1, mxREAL);
        }
        
        // Calculate the bispectrm / its gradient
        r_bispectrum(&shc1, r_lookup, &cgs, mxGetDoubles(plhs[0]));
        if (nlhs>1)
        {
            r_bispectrum_gradient(&shc1, r_lookup, &cgs, mxGetDoubles(plhs[1]));
        }
    }
    else
    {
        plhs[0] = mxCreateDoubleMatrix(c_lookup[bandlimit][bandlimit][bandlimit][1]+1, 1, mxREAL);
        shc2 = c_init_shc(bandlimit, mxGetDoubles(prhs[0]));

        if (nlhs>1)
        {
            plhs[1] = mxCreateDoubleMatrix(2*(bandlimit+1)*(bandlimit+1), c_lookup[bandlimit][bandlimit][bandlimit][1]+1, mxREAL);
        }
        
        // Calculate the bispectrm / its gradient
        c_bispectrum(&shc2, c_lookup, &cgs, mxGetDoubles(plhs[0]));
        if (nlhs>1)
        {
            c_bispectrum_gradient(&shc2, c_lookup, &cgs, mxGetDoubles(plhs[1]));
        }
    }
}
