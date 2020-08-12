/** 
 * Calculate the bispectrum and its gradient.
 * 
 * MATLAB call form:
 *      b = bispectrum_mex(shc, bandlimit)
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
// mxComplexDouble *shc_conj;
size_t bandlimit;
CGTable cgt;
mxArray *CGs;

/* Output variables */
double *b;
double *grad_i;
double *grad_j;
double *grad_vals;
long grad_nnz;
long grad_rows_no;
long grad_cols_no;


/* Utility variables */
blt lookup;
size_t grad_nnz; // Number of non-zero elements in the bispectrum gradient matrix
bool first_run = true;
size_t previous_bandlimit;


/* utility functions */
void build_bisp_lookup_table()
{
    grad_nnz = 0;
    
    size_t index = 0;

    lookup = malloc((bandlimit+1)*sizeof(size_t ***));

    for (long l1 = 0; l1<=bandlimit; ++l1)
    {
        lookup[l1] = malloc((l1+1)*sizeof(size_t **));
        for (long l2=0; l2<=l1; ++l2)
        {
            lookup[l1][l2] = malloc((((l1+l2)>=bandlimit ? bandlimit : l1+l2)  - (l1>=l2 ? l1-l2 : l2-l1) + 1)*sizeof(size_t *));
            for (long l = (l1>=l2 ? l1-l2 : l2-l1); l<=l1+l2 && l<=bandlimit; ++l)
            {
                lookup[l1][l2][l - (l1>=l2 ? l1-l2 : l2-l1)] = malloc(2*sizeof(size_t));
                lookup[l1][l2][l - (l1>=l2 ? l1-l2 : l2-l1)][0] = index++;
                lookup[l1][l2][l - (l1>=l2 ? l1-l2 : l2-l1)][1] = index++;

                if (l1==l2 && l2==l)
                {
                    grad_nnz += 1;
                }
                else if (l1==l2)
                {
                    grad_nnz += 2;
                }
                else if (l1==l)
                {
                    grad_nnz += 2;
                }
                else if (l2==l)
                {
                    grad_nnz += 2;
                }
                else
                {
                    grad_nnz += 3;
                }
            }
        }
    }

    grad_nnz *= 2;
    grad_rows_no = index;
    grad_cols_no = 2*(bandlimit+1)*(bandlimit+1);
}


void destory_bisp_lookup_table()
{
    for (long l1 = 0; l1<=bandlimit; ++l1)
    {
        for (long l2=0; l2<=l1; ++l2)
        {
            for (long l = (l1>=l2 ? l1-l2 : l2-l1); l<=l1+l2 && l<=bandlimit; ++l)
            {
                free(lookup[l1][l2][l - (l1>=l2 ? l1-l2 : l2-l1)]);
            }
            free(lookup[l1][l2]);
        }
        free(lookup[l1]);
    }
    free(lookup);
}


void create_CGTable()
{
    mxArray *temp[3];

    cgt = malloc((bandlimit+1)*sizeof(double ****));

    for (long l1 = 0; l1<=bandlimit; ++l1)
    {
        temp[0] = mxGetCell(CGs, l1);
        cgt[l1] = malloc((l1+1)*sizeof(double ***));

        for (long l2 = 0; l2<=l1; ++l2)
        {
            temp[1] = mxGetCell(temp[0], l2);
            cgt[l1][l2] = malloc(((l1+l2>=bandlimit ? bandlimit : l1+l2) - (l1-l2) + 1)*sizeof(double **));

            for (long l = l1-l2; l<=(l1+l2>=bandlimit ? bandlimit : l1+l2); ++l)
            {
                temp[2] = mxGetCell(temp[1], l-(l1-l2));
                cgt[l1][l2][l-(l1-l2)] = malloc((2*l+1)*sizeof(double *));

                for (long m = -l; m<=l; ++m)
                {
                    cgt[l1][l2][l-(l1-l2)][m+l] = mxGetDoubles(mxGetCell(temp[2], m+l));
                }
            }
        }
    }
}

void destroy_CGTable()
{
    for (long l1 = 0; l1<=bandlimit; ++l1)
    {
        for (long l2 = 0; l2<=l1; ++l2)
        {
            for (long l = l1-l2; l<=(l1+l2>=bandlimit ? bandlimit : l1+l2); ++l)
            {
                free(cgt[l1][l2][l-(l1-l2)]);
            }
            free(cgt[l1][l2]);
        }
        free(cgt[l1]);
    }
    free(cgt);
}

double get_cg(const long l1, const long l2, const long l, const long m, const long m1)
{
    cgt[l1][l2][l-(l1>=l2 ? l1-l2 : l2-l1)][m+l][m1-((m-l2<=-l1) ? -l1 : m-l2)];
}

mxComplexDouble get_shc(const long l, const long m)
{
    return shc[l*(l+1)+m];
}

long bisp_lookup(const long l1, const long l2, const long l, const PART part)
{
    return lookup[l1][l2][l - (l1-l2)][part];
}


/* Functions calculating the bispectrum and its gradient */

void bisp()
{
    for (long l1 = 0; l1<=bandlimit; ++l1)
    {
        for (long l2 = 0; l2<=l1; ++l2)
        {
            for (long l = l1-l2; l<=(l1+l2>=bandlimit ? bandlimit : l1+l2); ++l)
            {
                
                // TOOD
            }
        }
    }
}

void bisp_grad()
{
    // TODO
}


/**
 * The MEX Function
  */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Obtain input
    shc = mxGetComplexDoubles(prhs[0]);
    bandlimit = mxGetScalar(prhs[1]);

    // Handle first run
    if (first_run==true)
    {
        build_bisp_lookup_table();

        mxArray *inputMxArray[1];
        inputMxArray[0] = mxCreateDoubleScalar(bandlimit);
        mexCallMATLAB(0, NULL, 1, inputMxArray, "loadCGTable");

        CGs = mexGetVariablePtr("global", "CGs");
        create_CGTable();

        previous_bandlimit = bandlimit;
        first_run = false;
    }

    // Handle change of bandlimit
    if (bandlimit!=previous_bandlimit)
    {
        previous_bandlimit = bandlimit;

        destory_bisp_lookup_table();
        destroy_CGTable();

        mxArray *inputMxArray[1];
        inputMxArray[0] = mxCreateDoubleScalar(bandlimit);
        mexCallMATLAB(0, NULL, 1, inputMxArray, "loadCGTable");

        CGs = mexGetVariablePtr("global", "CGs");

        build_bisp_lookup_table();
        create_CGTable();
    }

    // Create output
    plhs[0] = mxCreateDoubleMatrix(grad_rows_no, 1, mxREAL);
    b = mxGetDoubles(plhs[0]);

    if (nlhs>1)
    {
        plhs[1] = mxCreateDoubleMatrix(grad_nnz, 1, mxREAL);
        grad_i = mxGetDoubles(plhs[1]);

        plhs[2] = mxCreateDoubleMatrix(grad_nnz, 1, mxREAL);
        grad_j = mxGetDoubles(plhs[2]);

        plhs[3] = mxCreateDoubleMatrix(grad_nnz, 1, mxREAL);
        grad_vals = mxGetDoubles(plhs[3]);

        plhs[4] = mxCreateDoubleScalar((double) grad_nnz);
        plhs[5] = mxCreateDoubleScalar((double) grad_rows_no);
        plhs[6] = mxCreateDoubleScalar((double) grad_cols_no);
    }

    // Compute output
    bisp();

    if (nlhs>1)
    {
        grad_bisp();
    }
}
