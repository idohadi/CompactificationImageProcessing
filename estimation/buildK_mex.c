/** 
 * TODO: 
 *      (1) Update it to use the new code base for Clbesch-Gordan coeffs
 *      (2) Write better docs
 * 
 * MATLAB call form:
 *      K1 = buildK_mex(UUH, UUT, bandlimit)
 *  UUH = U*(U^*) and UUT = U*(U^T)
 * 
 * NOTE:
 *  Code performs no input checks.
  */


#include <stdint.h>
#include "mex.h"
#include "clebsch_gordan_coefficients.h"

typedef enum PART {REAL_PART, IMAG_PART} PART;
typedef size_t **** c_bispectrum_lookout_table;
typedef c_bispectrum_lookout_table c_blt;

bool first_run = true;

c_blt lookup;
cg_table cgs;

size_t bandlimit = 0;
size_t previous_bandlimit = 0;

mxComplexDouble *UUH;
mxComplexDouble *UUT;

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
        lookup = c_build_bispectrum_lookup_table(bandlimit);
        mexPrintf("Initialization completed.\n");
        first_run = false;
    }

    // Handle chagne in bandlimit between runs
    if (previous_bandlimit!=bandlimit)
    {
        mexPrintf("Bandlimit chagned. Reinitializing Clebsch-Gordan coefficients and bispectrum lookup table.\n");
        // Destroying the lookup tables and Clebsch-Gordan table
        destroy_cg_table(&cgs);
        c_destroy_bispectrum_lookup_table(lookup, previous_bandlimit);

        // Generate new tables
        previous_bandlimit = bandlimit;
        cgs = allocate_cg_table(bandlimit);

        calculate_cg_table(bandlimit, &cgs);
        lookup = c_build_bispectrum_lookup_table(bandlimit);

        mexPrintf("Reinitialization completed.\n");
    }

    // TODO: actual code
}