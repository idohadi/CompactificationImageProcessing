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

mxDouble *K;

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
    UUH = mxGetComplexDoubles(prhs[0]);
    UUT = mxGetComplexDoubles(prhs[1]);
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

    // Generate output
    // TODO: improve sz,sz2 convention and use below, so as to make code clearer
    long sz = 2*(bandlimit+1)*(bandlimit+1);
    long sz2 = (bandlimit+1)*(bandlimit+1);
    plhs[0] = mxCreateDoubleMatrix(sz, lookup[bandlimit][bandlimit][bandlimit][1]+1, mxREAL);
    K = mxGetDoubles(plhs[0]);

    // Build K
    double rp, ip; // real part and imaginary part
    for (long l1 = 0; l1<=bandlimit; ++l1)
    {
        for (long l2 = 0; l2<=l1; ++l2)
        {
            for (long l = l1-l2; l<=bandlimit && l<=l1+l2; ++l)
            {
                // K1
                for (long m1 = -l1; m1<=l1; ++m1)
                {
                    rp = 0;
                    ip = 0;

                    for (long m = (-l>=m1-l2) ? -l : m1-l2; m<=l && m<=m1+l2; ++m)
                    {
                        rp += get_cg(&cgs, l1, l2, l, m, m1)*UUH[sz2*(l2*(l2+1) + (m-m1)) + (l*(l+1) + m)].real;
                        ip += get_cg(&cgs, l1, l2, l, m, m1)*UUH[sz2*(l2*(l2+1) + (m-m1)) + (l*(l+1) + m)].imag;
                    }

                    // Real part bisp invariant
                    K[sz*lookup[l1][l2][l - (l1-l2)][REAL_PART] + 2*(l1*(l1+1) + m1)] += rp;    // Real part SHC
                    K[sz*lookup[l1][l2][l - (l1-l2)][REAL_PART] + 2*(l1*(l1+1) + m1)+1] += ip;  // Imag part SHC
                    
                    // Imag part bisp invariant
                    K[sz*lookup[l1][l2][l - (l1-l2)][IMAG_PART] + 2*(l1*(l1+1) + m1)] += ip;    // Real part SHC
                    K[sz*lookup[l1][l2][l - (l1-l2)][IMAG_PART] + 2*(l1*(l1+1) + m1)+1] -= rp;   // Imag part SHC
                }

                // K2
                for (long m2 = -l2; m2<=l2; ++m2)
                {
                    rp = 0;
                    ip = 0;

                    for (long m = (-l>=m2-l1) ? -l : m2-l1; m<=l && m<=m2+l1; ++m)
                    {
                        rp += get_cg(&cgs, l1, l2, l, m, m-m2)*UUH[sz2*(l1*(l1+1) + (m-m2)) + (l*(l+1) + m)].real;
                        ip += get_cg(&cgs, l1, l2, l, m, m-m2)*UUH[sz2*(l1*(l1+1) + (m-m2)) + (l*(l+1) + m)].imag;
                    }

                    // Real part bisp invariant
                    K[sz*lookup[l1][l2][l - (l1-l2)][REAL_PART] + 2*(l2*(l2+1) + m2)] += rp;    // Real part SHC
                    K[sz*lookup[l1][l2][l - (l1-l2)][REAL_PART] + 2*(l2*(l2+1) + m2)+1] += ip;  // Imag part SHC
                    
                    // Imag part bisp invariant
                    K[sz*lookup[l1][l2][l - (l1-l2)][IMAG_PART] + 2*(l2*(l2+1) + m2)] += ip;    // Real part SHC
                    K[sz*lookup[l1][l2][l - (l1-l2)][IMAG_PART] + 2*(l2*(l2+1) + m2)+1] -= rp;   // Imag part SHC
                }

                // K3
                for (long m = -l; m<=l; ++m)
                {
                    rp = 0;
                    ip = 0;

                    for (long m1 = (-l1>=m-l2) ? -l1 : m-l2; m1<=l1 && m1<=m+l2; ++m1)
                    {
                        rp += get_cg(&cgs, l1, l2, l, m, m1)*UUT[sz2*(l2*(l2+1) + (m-m1)) + (l1*(l1+1) + m1)].real;
                        ip -= get_cg(&cgs, l1, l2, l, m, m1)*UUT[sz2*(l2*(l2+1) + (m-m1)) + (l1*(l1+1) + m1)].imag;
                    }

                    // Real part bisp invariant
                    K[sz*lookup[l1][l2][l - (l1-l2)][REAL_PART] + 2*(l*(l+1) + m)] += rp;    // Real part SHC
                    K[sz*lookup[l1][l2][l - (l1-l2)][REAL_PART] + 2*(l*(l+1) + m)+1] -= ip;  // Imag part SHC
                    
                    // Imag part bisp invariant
                    K[sz*lookup[l1][l2][l - (l1-l2)][IMAG_PART] + 2*(l*(l+1) + m)] += ip;    // Real part SHC
                    K[sz*lookup[l1][l2][l - (l1-l2)][IMAG_PART] + 2*(l*(l+1) + m)+1] += rp;   // Imag part SHC
                }

            }
        }
    }

}