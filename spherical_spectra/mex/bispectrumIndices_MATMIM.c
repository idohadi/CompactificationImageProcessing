// TODO: docs in doxygen
/** 
 * Return a bispectrum vector based on spherical harmonics vector. 
 * 
 * MATLAB call form:
 *      inds = bispectrumIndices_MATMIM(bandlimit, l_target)
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

typedef enum PART {REAL_PART, IMAG_PART} PART;
typedef size_t **** c_bispectrum_lookout_table;
typedef c_bispectrum_lookout_table c_blt;

bool first_run = true;

c_blt lookup;

size_t bandlimit = 0;
size_t previous_bandlimit = 0;
size_t l_target = 0;

double *output;

c_blt c_build_bispectrum_lookup_table(size_t bandlimit)
{
    size_t index = 0; 

    size_t ****lookup_table = malloc((bandlimit+1)*sizeof(size_t ***));
    for (long l1 = 0; l1<=bandlimit; ++l1)
    {
        lookup_table[l1] = malloc((l1+1)*sizeof(size_t **));
        for (long l2=0; l2<=l1; ++l2)
        {
            lookup_table[l1][l2] = malloc((((l1+l2)>=bandlimit ? bandlimit : l1+l2)  - (l1>=l2 ? l1-l2 : l2-l1) + 1)*sizeof(size_t *));
            for (long l = (l1>=l2 ? l1-l2 : l2-l1); l<=l1+l2 && l<=bandlimit; ++l)
            {
                lookup_table[l1][l2][l - (l1>=l2 ? l1-l2 : l2-l1)] = malloc(2*sizeof(size_t));
                lookup_table[l1][l2][l - (l1>=l2 ? l1-l2 : l2-l1)][0] = index++;
                lookup_table[l1][l2][l - (l1>=l2 ? l1-l2 : l2-l1)][1] = index++;
            }
        }
    }
    return lookup_table;
}

void c_destroy_bispectrum_lookup_table(c_blt table, size_t bandlimit)
{
    for (long l1 = 0; l1<=bandlimit; ++l1)
    {
        for (long l2=0; l2<=l1; ++l2)
        {
            for (long l = (l1>=l2 ? l1-l2 : l2-l1); l<=l1+l2 && l<=bandlimit; ++l)
            {
                free(table[l1][l2][l - (l1>=l2 ? l1-l2 : l2-l1)]);
            }
            free(table[l1][l2]);
        }
        free(table[l1]);
    }
    free(table);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Get input variables
    bandlimit = (size_t) mxGetScalar(prhs[0]);
    l_target = (size_t) mxGetScalar(prhs[1]);

    // Handle first run
    if (first_run == true)
    {
        mexPrintf("First run. Initializing bispectrum lookup table.\n");
        previous_bandlimit = bandlimit;
        lookup = c_build_bispectrum_lookup_table(bandlimit);
        mexPrintf("Initialization completed.\n");
        first_run = false;
    }

    // Handle chagne in bandlimit between runs
    if (previous_bandlimit!=bandlimit)
    {
        mexPrintf("Bandlimit chagned. Reinitializing bispectrum lookup table.\n");
        // Destroying the lookup tables and Clebsch-Gordan table
        c_destroy_bispectrum_lookup_table(lookup, previous_bandlimit);

        // Generate new table
        previous_bandlimit = bandlimit;
        lookup = c_build_bispectrum_lookup_table(bandlimit);

        mexPrintf("Reinitialization completed.\n");
    }

    // Create output
    long count = 0;
    for (long l1 = 0; l1<=bandlimit && l1<=l_target; ++l1)
    {
        for (long l2 = 0; l2<=l1 && l2<=l_target; ++l2)
        {
            for (long l = l1-l2; l<=bandlimit && l<=l1+l2 && l<=l_target; ++l)
            {
                if (l1 == l_target)
                {
                    ++count;
                } 
                else if (l2 == l_target)
                {
                    ++count;
                }
                else if (l == l_target)
                {
                    ++count;
                }
            }
        }
    }

    plhs[0] = mxCreateDoubleMatrix(2*count, 1, mxREAL);
    output = mxGetDoubles(plhs[0]);

    count = 0;
    for (long l1 = 0; l1<=bandlimit && l1<=l_target; ++l1)
    {
        for (long l2 = 0; l2<=l1 && l2<=l_target; ++l2)
        {
            for (long l = l1-l2; l<=bandlimit && l<=l1+l2 && l<=l_target; ++l)
            {
                if (l1 == l_target)
                {
                    output[count++] = lookup[l1][l2][l-(l1-l2)][REAL_PART] + 1;
                    output[count++] = lookup[l1][l2][l-(l1-l2)][IMAG_PART] + 1;
                } 
                else if (l2 == l_target)
                {
                    output[count++] = lookup[l1][l2][l-(l1-l2)][REAL_PART] + 1;
                    output[count++] = lookup[l1][l2][l-(l1-l2)][IMAG_PART] + 1;
                }
                else if (l == l_target)
                {
                    output[count++] = lookup[l1][l2][l-(l1-l2)][REAL_PART] + 1;
                    output[count++] = lookup[l1][l2][l-(l1-l2)][IMAG_PART] + 1;
                }
            }
        }
    }
}
