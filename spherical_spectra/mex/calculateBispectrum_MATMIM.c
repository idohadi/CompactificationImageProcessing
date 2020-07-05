// TODO: docs in doxygen
/** 
 * Return a bispectrum vector based on spherical harmonics vector. 
 * 
 * MATLAB call form:
 *      [bispectrum, bispectrum_gradient] = calculateBispectrum_MATMIM(shc, shc_conjugated, bandlimit)
 *  where
 *      bispectrum is the bispectrum vector. 
 *      bispectrum_gradient (optional argument) is its gradient (its transpose, actually).
 *      shc is assumed to represent complex-valued spherical function. It is a complex row or column vector.
 *      shc_conjugated is the complex conjugate of shc.
 * 
 * NOTE:
 *  Code performs no input checks.
  */


#include <stdint.h>
#include "mex.h"
#include "spherical_spectra.h"

bool first_run = true;
c_blt c_lookup;
cg_table cgs;
size_t bandlimit = 0;
size_t previous_bandlimit = 0;
mxComplexDouble *shc;
mxComplexDouble *shc_conjugated;

void complexSum(mxComplexDouble *arg1, mxComplexDouble *arg2, mxComplexDouble *output)
{
    output->real = arg1->real + arg2->real;
    output->imag = arg1->imag + arg2->imag;
}

void complexProduct(mxComplexDouble *arg1, mxComplexDouble *arg2, mxComplexDouble *output)
{
    output->real = arg1->real*arg2->real - arg1->imag*arg2->imag;
    output->imag = arg1->real*arg2->imag + arg1->imag*arg2->real;
}

void bispectral_invariant(  mxComplexDouble * const shc, mxComplexDouble * const shc_conjugated, 
                            const long l1, const long l2, const long l, const size_t bandlimit, 
                            const c_blt lookup, cg_table * const table, 
                            double *output_rp, double *output_ip)
{
    // TODO: calculate the bispectral invariant l1, l2, l
}

void bisp(  mxComplexDouble * const shc, mxComplexDouble * const shc_conjugated, const size_t bandlimit, 
            const c_blt lookup, cg_table * const table, double *output)
{
    // TODO: calculate all the bispectral invariants
}

void bispGrad(  mxComplexDouble * const shc, mxComplexDouble * const shc_conjugated, const size_t bandlimit, 
                const c_blt lookup, cg_table * const table, double *output)
{
    // TODO: calculate all the gradients of bispectral invariants 
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Get input variables
    shc = mxGetComplexDoubles(prhs[0]);
    shc_conjugated = mxGetComplexDoubles(prhs[1]);
    bandlimit = (size_t) mxGetScalar(prhs[2]);

    // Handle first run
    if (first_run == true)
    {
        mexPrintf("First run. Initializing Clebsch-Gordan coefficients and bispectrum lookup table.\n");
        previous_bandlimit = bandlimit;
        cgs = allocate_cg_table(bandlimit);
        calculate_cg_table(bandlimit, &cgs);
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
        c_destroy_bispectrum_lookup_table(c_lookup, previous_bandlimit);

        // Generate new tables
        previous_bandlimit = bandlimit;
        cgs = allocate_cg_table(bandlimit);

        calculate_cg_table(bandlimit, &cgs);
        c_lookup = c_build_bispectrum_lookup_table(bandlimit);

        mexPrintf("Reinitialization completed.\n");
    }

    // Create output
    plhs[0] = mxCreateDoubleMatrix(c_lookup[bandlimit][bandlimit][bandlimit][1]+1, 1, mxREAL);
    if (nlhs>1)
    {
        plhs[1] = mxCreateDoubleMatrix(2*(bandlimit+1)*(bandlimit+1), c_lookup[bandlimit][bandlimit][bandlimit][1]+1, mxREAL);
    }
    
    // Calculate the bispectrm / its gradient
    bisp(shc, shc_conjugated, bandlimit, c_lookup, &cgs, mxGetDoubles(plhs[0]));
    if (nlhs>1)
    {
        bispGrad(shc, shc_conjugated, bandlimit, c_lookup, &cgs, mxGetDoubles(plhs[1]));
    }
}