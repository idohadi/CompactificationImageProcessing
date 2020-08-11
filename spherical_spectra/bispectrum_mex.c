/** 
 * Calculate the bispectrum and its gradient.
 * 
 * MATLAB call form:
 *      b = bispectrum_mex(shc, shc_conj, bandlimit, CGs)
 *      [b, grad_i, grad_j, grad_vals, grad_nnz, grad_rows_no, grad_cols_no] = bispectrum_mex(shc, shc_conj, bandlimit, CGs)
 *  where
 *      shc                           column or row complex array of length 
 *                                    (bandlimit+1)^2 of spherical harmonics coefficients
 *      shc_conj                      the complex conjugate of shc
 *      bandlimit                     scalar, the bandlimit of the function represented 
 *                                    by shc
 *      CGs                           a cell array containing the precomputed Clebsch-Gordan
 *                                    coefficents for bandlimit. 
 *                                    Their format is documented bispectrum.m docs.
 *      b                             bispectrum column vector, containing the bispectrum 
 *                                    invariants b_{l1,l2,l} for 
 *                                         0<=l1<=bandlimit, 0<=l2<=l1, abs(l1-l2)<=l<=min(bandlmit, l1+l2)
 *                                    separated to real and complex parts and ordered in lexigraphical order.
 *                                    That is, 
 *                                      b = [real(b_{0,0,0}); imag(b_{0,0,0}); ...
 *                                          real(b_{1,0,0}); imag(b_{1,0,0}); ...]
 *      grad_i, grad_j, grad_vals      grad_nnz x 1 arrays, such that the gardient of the bispectrum satisfies
 *                                        grad(grad_i(k), grad_j(k)) = grad_vals(k)
 *      grad_nnz                      number of non-zero elements in the gradient
 *      grad_rows_no, grad_cols_no    number of rows and columns of the gradient
 * 
 * 
 * NOTES:
 *  (1) The code performs no input checks.
 *  (2) This function is wrapped by bispectrum.m.
 * 
  */

#include <stdint.h>
#include "mex.h"


/* Typedefs */
typedef enum PART {REAL_PART, IMAG_PART} PART;

typedef size_t **** bispectrum_lookup_table;
typedef bispectrum_lookup_table blt;

typedef double ***** CGTable;


/* Input variables */
mxComplexDouble *shc;
mxComplexDouble *shc_conj;
size_t bandlimit;
CGTable cg;


/* Output variables */
double *b;
double *grad_inds;
double *grad_vals;


/* Utility variables */
blt lookup;
size_t grad_nnz; // Number of non-zero elements in the bispectrum gradient matrix


/* utility functions */

void ccsum(mxComplexDouble *arg1, mxComplexDouble *arg2, mxComplexDouble *output)
{
  /** 
   * Perform output = arg1 + arg2 for complex variables.
   * */

    output->real = arg1->real + arg2->real;
    output->imag = arg1->imag + arg2->imag;
}

void ccprod_addinsitu(mxComplexDouble *arg1, mxComplexDouble *arg2, double factor, mxComplexDouble *output)
{
    /** 
   * Perform output += factor*arg1*arg2 for complex variables.
   * */

    output->real += factor*(arg1->real*arg2->real - arg1->imag*arg2->imag);
    output->imag += factor*(arg1->real*arg2->imag + arg1->imag*arg2->real);
}

void rcprod(double *arg1, mxComplexDouble *arg2, mxComplexDouble *output)
{
  /**
   * Perform output = arg1*arg2 for complex arg2 and real arg1.
   */

    output->real = (*arg1)*arg2->real;
    output->imag = (*arg1)*arg2->imag;
}

void build_bisp_lookup_table()
{
  // TODO
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


void destory_bisp_lookup_table(c_blt table, size_t bandlimit)
{
  // TODO
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


void create_CGTable()
{
  // TODO
}

void destroy_CGTable()
{
  // TODO
}


/* Functions calculating the bispectrum and its gradient */

void bisp(mxComplexDouble * const input, mxComplexDouble * const input_conjugated, double *output)
{
  // TODO
    mxComplexDouble temp;

    for (long l1 = 0; l1<=bandlimit; ++l1)
    {
        for (long l2 = 0; l2<=l1; ++l2)
        {
            for (long l = l1-l2; l<=bandlimit && l<=l1+l2; ++l)
            {
                bispectral_invariant(input, input_conjugated, l1, l2, l, &temp);
                output[lookup[l1][l2][l - (l1-l2)][REAL_PART]] = temp.real;
                output[lookup[l1][l2][l - (l1-l2)][IMAG_PART]] = temp.imag;
            }
        }
    }
}

void bispGrad(  mxComplexDouble * const input, mxComplexDouble * const input_conjugated, double *output)
{
    // TODO
    
    mxArray *arr = mxCreateDoubleMatrix((bandlimit+1)*(bandlimit+1), 1, mxCOMPLEX);
    mxComplexDouble *onei = mxGetComplexDoubles(arr);
    
    mxComplexDouble temp;

    size_t row_index_rp, row_index_ip;
    
    for (long l1 = 0; l1<=bandlimit; ++l1)
    {
        for (long l2 = 0; l2<=l1; ++l2)
        {
            for (long l = l1-l2; l<=bandlimit && l<=l1+l2; ++l)
            {
                row_index_rp = 2*(bandlimit+1)*(bandlimit+1)*lookup[l1][l2][l-(l1-l2)][REAL_PART];
                row_index_ip = 2*(bandlimit+1)*(bandlimit+1)*lookup[l1][l2][l-(l1-l2)][IMAG_PART];

                // Derivative of b_{l1,l2,l} w.r.t f_{l1,m1_der}
                for (long m1_der = -l1; m1_der<=l1; ++m1_der)
                {
                    // Derivative w.r.t real part
                    temp.real = 0;
                    temp.imag = 0;
                    onei[l1*(l1+1) + m1_der].real = 1;
                    sum_defining_bisp(onei, shc_conjugated, shc, l1, l2, l, &temp);
                    output[row_index_rp + 2*(l1*(l1+1) + m1_der)] += temp.real; // Derivative of real part of b_{l1,l2,l}
                    output[row_index_ip + 2*(l1*(l1+1) + m1_der)] += temp.imag; // Derivative of imag part of b_{l1,l2,l}
                    onei[l1*(l1+1) + m1_der].real = 0;

                    // Derivative w.r.t imaginary part
                    temp.real = 0;
                    temp.imag = 0;
                    onei[l1*(l1+1) + m1_der].imag = -1;
                    sum_defining_bisp(onei, shc_conjugated, shc, l1, l2, l, &temp);
                    output[row_index_rp + 2*(l1*(l1+1) + m1_der)+1] += temp.real; // Derivative of real part of b_{l1,l2,l}
                    output[row_index_ip + 2*(l1*(l1+1) + m1_der)+1] += temp.imag; // Derivative of imag part of b_{l1,l2,l}
                    onei[l1*(l1+1) + m1_der].imag = 0;
                }


                // Derivative of b_{l1,l2,l} w.r.t f_{l2,m2_der}
                for (long m2_der = -l2; m2_der<=l2; ++m2_der)
                {
                    // Derivative w.r.t real part
                    temp.real = 0;
                    temp.imag = 0;
                    onei[l2*(l2+1) + m2_der].real = 1;
                    sum_defining_bisp(shc_conjugated, onei, shc, l1, l2, l, &temp);
                    output[row_index_rp + 2*(l2*(l2+1) + m2_der)] += temp.real; // Derivative of real part of b_{l1,l2,l}
                    output[row_index_ip + 2*(l2*(l2+1) + m2_der)] += temp.imag; // Derivative of imag part of b_{l1,l2,l}
                    onei[l2*(l2+1) + m2_der].real = 0;

                    // Derivative w.r.t imaginary part
                    temp.real = 0;
                    temp.imag = 0;
                    onei[l2*(l2+1) + m2_der].imag = -1;
                    sum_defining_bisp(shc_conjugated, onei, shc, l1, l2, l, &temp);
                    output[row_index_rp + 2*(l2*(l2+1) + m2_der)+1] += temp.real; // Derivative of real part of b_{l1,l2,l}
                    output[row_index_ip + 2*(l2*(l2+1) + m2_der)+1] += temp.imag; // Derivative of imag part of b_{l1,l2,l}
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
                    output[row_index_rp + 2*(l*(l+1) + m_der)] += temp.real; // Derivative of real part of b_{l1,l2,l}
                    output[row_index_ip + 2*(l*(l+1) + m_der)] += temp.imag; // Derivative of imag part of b_{l1,l2,l}
                    onei[l*(l+1) + m_der].real = 0;

                    // Derivative w.r.t imaginary part
                    temp.real = 0;
                    temp.imag = 0;
                    onei[l*(l+1) + m_der].imag = 1;
                    sum_defining_bisp(shc_conjugated, shc_conjugated, onei, l1, l2, l, &temp);
                    output[row_index_rp + 2*(l*(l+1) + m_der)+1] += temp.real; // Derivative of real part of b_{l1,l2,l}
                    output[row_index_ip + 2*(l*(l+1) + m_der)+1] += temp.imag; // Derivative of imag part of b_{l1,l2,l}
                    onei[l*(l+1) + m_der].imag = 0;
                }
            }
        }
    }
}

/**
 * The MEX Function
  */
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
