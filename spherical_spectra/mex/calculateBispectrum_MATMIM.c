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
#include "clebsch_gordan_coefficients.h"
// #include "spherical_spectra.h"

typedef enum PART {REAL_PART, IMAG_PART} PART;
typedef size_t **** c_bispectrum_lookout_table;
typedef c_bispectrum_lookout_table c_blt;

bool first_run = true;

c_blt lookup;
cg_table cgs;

size_t bandlimit = 0;
size_t previous_bandlimit = 0;

mxComplexDouble *shc;
mxComplexDouble *shc_conjugated;

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

void ccsum(mxComplexDouble *arg1, mxComplexDouble *arg2, mxComplexDouble *output)
{
    output->real = arg1->real + arg2->real;
    output->imag = arg1->imag + arg2->imag;
}

void ccprod_addinsitu(mxComplexDouble *arg1, mxComplexDouble *arg2, double factor, mxComplexDouble *output)
{
    // Perform output += factor*arg1*arg2

    output->real += factor*(arg1->real*arg2->real - arg1->imag*arg2->imag);
    output->imag += factor*(arg1->real*arg2->imag + arg1->imag*arg2->real);
}

void rcprod(double *arg1, mxComplexDouble *arg2, mxComplexDouble *output)
{
    output->real = (*arg1)*arg2->real;
    output->imag = (*arg1)*arg2->imag;
}

mxComplexDouble *get_shc_ptr(mxComplexDouble * const input, const long l, const long m)
{
    if (l > bandlimit)
   {
       return NULL;
   }
   else if (m<-l || m>l)
   {
       mexPrintf("SHC (l,m) indices are non-standard.\n");
       return NULL;
   }
   else
   {
       return &input[l*(l + 1) + m];
   }
}

void bispectral_invariant(  mxComplexDouble * const input, mxComplexDouble * const input_conjugated, 
                            const long l1, const long l2, const long l, 
                            mxComplexDouble *output)
{
    output->real = 0;
    output->imag = 0;

    mxComplexDouble temp;

    for (long m = - l; m<=l; ++m)
    {
        temp.real = 0;
        temp.imag = 0;
        for (long m1 = ((m-l2)<=-l1 ? -l1 : (m-l2)); m1 <= ((m+l2)<=l1 ? (m+l2) : l1); ++m1)
        {
            ccprod_addinsitu(   get_shc_ptr(input_conjugated, l1, m1), 
                                get_shc_ptr(input_conjugated, l2, m-m1), 
                                get_cg(&cgs, l1, l2, l, m, m1), 
                                &temp);
        }
        ccprod_addinsitu(get_shc_ptr(input, l, m), &temp, 1.0, output);
    }
}

// mxComplexDouble conj(mxComplexDouble * const arg)
// {
//     mxComplexDouble output;
//     output.real = arg->real;
//     output.imag = -arg->imag;
//     return output;
// }

void bisp(mxComplexDouble * const input, mxComplexDouble * const input_conjugated, double *output)
{
    mxComplexDouble temp;

    for (long l1 = 0; l1<=bandlimit; ++l1)
    {
        for (long l2 = 0; l2<=l1; ++l2)
        {
            for (long l = l1-l2; l<=((bandlimit>=l1+l2) ?  l1+l2 : bandlimit); ++l)
            {
                bispectral_invariant(input, input_conjugated, l1, l2, l, &temp);
                output[lookup[l1][l2][l - (l1-l2)][REAL_PART]] = temp.real;
                output[lookup[l1][l2][l - (l1-l2)][IMAG_PART]] = temp.imag;
            }
        }
    }
}

void sum_defining_bisp( mxComplexDouble * const input1, mxComplexDouble * const input2, mxComplexDouble * const input3, 
                        const long l1, const long l2, const long l, mxComplexDouble *output)
{
    // Calculate 
    //  sum_{|m1|<=l1, |m2|<=l2, |m|<=l} <l1 m1 l2 m2 | l m> * conj(input1_{l1,m1}) * conj(input2_{l2,m2}) * conj(input3_{l,m})
    // input1, input2, input3 are assumed to have size (bandlimit+1)^2
    // 
    output->real = 0;
    output->imag = 0;

    mxComplexDouble temp;

    for (long m = - l; m<=l; ++m)
    {
        temp.real = 0;
        temp.imag = 0;
        for (long m1 = ((m-l2)<=-l1 ? -l1 : (m-l2)); m1 <= ((m+l2)<=l1 ? (m+l2) : l1); ++m1)
        {
            ccprod_addinsitu(   get_shc_ptr(input1, l1, m1), 
                                get_shc_ptr(input2, l2, m-m1), 
                                get_cg(&cgs, l1, l2, l, m, m1), 
                                &temp);
        }
        ccprod_addinsitu(get_shc_ptr(input3, l, m), &temp, 1.0, output);
    }
}


void bispGrad(  mxComplexDouble * const input, mxComplexDouble * const input_conjugated, double *output)
{
    // TODO: calculate all the gradients of bispectral invariants 
    mxArray *arr = mxCreateDoubleMatrix((bandlimit+1)*(bandlimit+1), 1, mxCOMPLEX);
    mxComplexDouble *onei = mxGetComplexDoubles(arr);
    for (long i = 0; i<=(bandlimit+1)*(bandlimit+1); ++i)
    {
        onei[i].real = 1;
    }
    mxComplexDouble temp;

    size_t row_index_rp, row_index_ip;
    
    for (long l1 = 0; l1<=bandlimit; ++l1)
    {
        for (long l2 = 0; l2<=l1; ++l2)
        {
            for (long l = l1-l2; l<=((bandlimit>=l1+l2) ?  l1+l2 : bandlimit); ++l)
            {
                row_index_rp = 2*(bandlimit+1)*(bandlimit+1)*lookup[l1][l2][l-(l1-l2)][REAL_PART];
                row_index_ip = 2*(bandlimit+1)*(bandlimit+1)*lookup[l1][l2][l-(l1-l2)][IMAG_PART];

                // Derivative of b_{l1,l2,l} w.r.t f_{l1,m1_der}
                for (long m1_der = -l1; m1_der<=l1; ++m1_der)
                {
                    // Derivative w.r.t real part
                    temp.real = 0;
                    temp.imag = 0;
                    onei[l1*(l1+1) + m1_der].real = -1;
                    sum_defining_bisp(onei, shc_conjugated, shc, l1, l2, l, &temp);
                    output[row_index_rp + 2*(l1*(l1+1) + m1_der)] = temp.real; // Derivative of real part of b_{l1,l2,l}
                    output[row_index_ip + 2*(l1*(l1+1) + m1_der)] = temp.imag; // Derivative of imag part of b_{l1,l2,l}
                    onei[l1*(l1+1) + m1_der].real = 0;

                    // Derivative w.r.t imaginary part
                    temp.real = 0;
                    temp.imag = 0;
                    onei[l1*(l1+1) + m1_der].imag = -1;
                    sum_defining_bisp(onei, shc_conjugated, shc, l1, l2, l, &temp);
                    output[row_index_rp + 2*(l1*(l1+1) + m1_der)+1] = temp.real; // Derivative of real part of b_{l1,l2,l}
                    output[row_index_ip + 2*(l1*(l1+1) + m1_der)+1] = temp.imag; // Derivative of imag part of b_{l1,l2,l}
                    onei[l1*(l1+1) + m1_der].imag = 0;
                }


                // Derivative of b_{l1,l2,l} w.r.t f_{l2,m2_der}
                for (long m2_der = -l2; m2_der<=l2; ++m2_der)
                {
                    // Derivative w.r.t real part
                    temp.real = 0;
                    temp.imag = 0;
                    onei[l2*(l2+1) + m2_der].real = -1;
                    sum_defining_bisp(shc_conjugated, onei, shc, l1, l2, l, &temp);
                    output[row_index_rp + 2*(l2*(l2+1) + m2_der)] = temp.real; // Derivative of real part of b_{l1,l2,l}
                    output[row_index_ip + 2*(l2*(l2+1) + m2_der)] = temp.imag; // Derivative of imag part of b_{l1,l2,l}
                    onei[l2*(l2+1) + m2_der].real = 0;

                    // Derivative w.r.t imaginary part
                    temp.real = 0;
                    temp.imag = 0;
                    onei[l2*(l2+1) + m2_der].imag = -1;
                    sum_defining_bisp(shc_conjugated, onei, shc, l1, l2, l, &temp);
                    output[row_index_rp + 2*(l2*(l2+1) + m2_der)+1] = temp.real; // Derivative of real part of b_{l1,l2,l}
                    output[row_index_ip + 2*(l2*(l2+1) + m2_der)+1] = temp.imag; // Derivative of imag part of b_{l1,l2,l}
                    onei[l2*(l2+1) + m2_der].imag = 0;
                }


                // Derivative of b_{l1,l2,l} w.r.t f_{l,m_der}
                for (long m_der = -l; m_der<=l; ++m_der)
                {
                    // Derivative w.r.t real part
                    temp.real = 0;
                    temp.imag = 0;
                    onei[l*(l+1) + m_der].real = 1;
                    sum_defining_bisp(shc_conjugated, shc_conjugated, onei, l1, l2, l, &temp);
                    output[row_index_rp + 2*(l*(l+1) + m_der)] = temp.real; // Derivative of real part of b_{l1,l2,l}
                    output[row_index_ip + 2*(l*(l+1) + m_der)] = temp.imag; // Derivative of imag part of b_{l1,l2,l}
                    onei[l*(l+1) + m_der].real = 0;

                    // Derivative w.r.t imaginary part
                    temp.real = 0;
                    temp.imag = 0;
                    onei[l*(l+1) + m_der].imag = 1;
                    sum_defining_bisp(shc_conjugated, shc_conjugated, onei, l1, l2, l, &temp);
                    output[row_index_rp + 2*(l*(l+1) + m_der)+1] = temp.real; // Derivative of real part of b_{l1,l2,l}
                    output[row_index_ip + 2*(l*(l+1) + m_der)+1] = temp.imag; // Derivative of imag part of b_{l1,l2,l}
                    onei[l*(l+1) + m_der].imag = 0;
                }
            }
        }
    }
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

    // Create output
    plhs[0] = mxCreateDoubleMatrix(lookup[bandlimit][bandlimit][bandlimit][1]+1, 1, mxREAL);
    if (nlhs>1)
    {
        plhs[1] = mxCreateDoubleMatrix(2*(bandlimit+1)*(bandlimit+1), lookup[bandlimit][bandlimit][bandlimit][1]+1, mxREAL);
    }
    
    // Calculate the bispectrm / its gradient
    bisp(shc, shc_conjugated, mxGetDoubles(plhs[0]));
    if (nlhs>1)
    {
        bispGrad(shc, shc_conjugated, mxGetDoubles(plhs[1]));
    }
}