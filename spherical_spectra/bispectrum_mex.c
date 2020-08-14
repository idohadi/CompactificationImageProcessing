/** 
 * Calculate the bispectrum and its gradient.
 * 
 * MATLAB call form:
 *      b = bispectrum_mex(shc, bandlimit)
 *      [b, grad_i, grad_j, grad_vals, grad_nnz, grad_rows_no, grad_cols_no] = bispectrum_mex(shc, bandlimit)
 *  where
 *      shc                           column or row complex array of length 
 *                                    (bandlimit+1)^2 of spherical harmonics coefficients
 *      bandlimit                     scalar, the bandlimit of the function represented 
 *                                    by shc
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
const mxArray *CGs;

/* Output variables */
double *b;
double *grad_i;
double *grad_j;
double *grad_vals;
long grad_nnz; // Number of non-zero elements in the bispectrum gradient matrix
long grad_rows_no;
long grad_cols_no;


/* Utility variables */
blt lookup;
bool first_run = true;
size_t previous_bandlimit;

size_t nnz_ind;


/* utility functions */
void build_bisp_lookup_table()
{
    grad_nnz = 0;
    
    size_t index = 0;

    lookup = malloc((bandlimit+1)*sizeof(size_t ***));

    for (long l1 = 0; l1<=bandlimit; ++l1)
    {
        lookup[l1] = malloc((l1+1)*sizeof(size_t **));
        for (long l2 = 0; l2<=l1; ++l2)
        {
            lookup[l1][l2] = malloc((((l1+l2)>=bandlimit ? bandlimit : l1+l2)  - (l1-l2) + 1)*sizeof(size_t *));
            for (long l = l1-l2; l<=l1+l2 && l<=bandlimit; ++l)
            {
                lookup[l1][l2][l - (l1-l2)] = malloc(2*sizeof(size_t));
                lookup[l1][l2][l - (l1-l2)][0] = index++;
                lookup[l1][l2][l - (l1-l2)][1] = index++;

                if (l1==l2 && l2==l)
                {
                    grad_nnz += 2*l+1;
                }
                else if (l1==l2)
                {
                    grad_nnz += 2*(l1 + l) + 2;
                }
                else if (l1==l)
                {
                    grad_nnz += 2*(l1 + l2) + 2;
                }
                else if (l2==l)
                {
                    grad_nnz += 2*(l2 + l1) + 2;
                }
                else
                {
                    grad_nnz += 2*(l1+l2+l) + 3;
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
    for (long l1 = 0; l1<=previous_bandlimit; ++l1)
    {
        for (long l2 = 0; l2<=l1; ++l2)
        {
            for (long l = l1-l2; l<=l1+l2 && l<=previous_bandlimit; ++l)
            {
                free(lookup[l1][l2][l - (l1-l2)]);
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

            for (long l = l1-l2; l<=l1+l2 && l<=bandlimit; ++l)
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
    for (long l1 = 0; l1<=previous_bandlimit; ++l1)
    {
        for (long l2 = 0; l2<=l1; ++l2)
        {
            for (long l = l1-l2; l<=l1+l2 && l<=previous_bandlimit; ++l)
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
    return cgt[l1][l2][l-(l1>=l2 ? l1-l2 : l2-l1)][m+l][m1-((m-l2<=-l1) ? -l1 : m-l2)];
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
            for (long l = l1-l2; l<=l1+l2 && l<=bandlimit; ++l)
            {
                long ind_real = bisp_lookup(l1, l2, l, REAL_PART);
                long ind_imag = bisp_lookup(l1, l2, l, IMAG_PART);

                double b_real;
                double b_imag;
                for (long m = -l; m<=l; ++m)
                {
                    b_real = 0;
                    b_imag = 0;
                    
                    for (long m1 = ((m-l2)<=-l1 ? -l1 : (m-l2)); m1 <= ((m+l2)<=l1 ? (m+l2) : l1); ++m1)
                    {
                        b_real += get_cg(l1, l2, l, m, m1) 
                                    * (get_shc(l1, m1).real * get_shc(l2, m-m1).real 
                                        - get_shc(l1, m1).imag * get_shc(l2, m-m1).imag);
                        b_imag += - get_cg(l1, l2, l, m, m1) 
                                    * (get_shc(l1, m1).real * get_shc(l2, m-m1).imag 
                                        + get_shc(l1, m1).imag * get_shc(l2, m-m1).real);
                    }
                    
                    // b[ind_real] += get_shc(l, m).real * b_real - get_shc(l, m).imag * b_imag;
                    // b[ind_imag] += get_shc(l, m).imag * b_real + get_shc(l, m).real * b_imag;
                    b[ind_real] += b_real;
                    b[ind_imag] += b_imag;
                }
            }
        }
    }
}

long maximum(const long a, const long b)
{
    return (b>=a ? b : a);
}

long minimum(const long a, const long b)
{
    return (b>=a ? a : b);
}


void calc_l_sum(const long l1, const long l2, const long l, const long n, double *l_sum_real, double *l_sum_imag)
{
    for (long m1 = maximum(-l1, n-l2); m1<=minimum(l1, n+l2); ++m1)
    {
        *l_sum_real += get_cg(l1, l2, l, n, m1) *
                        (get_shc(l1, m1).real * get_shc(l2, n-m1).real 
                        - get_shc(l1, m1).imag * get_shc(l2, n-m1).imag);
        *l_sum_imag += -get_cg(l1, l2, l, n, m1) *
                        (get_shc(l1, m1).imag * get_shc(l2, n-m1).real
                        + get_shc(l1, m1).real * get_shc(l2, n-m1).imag);
    }
}


void calc_l1_sum(const long l1, const long l2, const long l, const long n, double *l1_sum_real, double *l1_sum_imag)
{
    for (long m = maximum(-l, n-l2); m<=minimum(l, n+l2); ++m)
    {
        *l1_sum_real += get_cg(l1, l2, l, m, n) *
                        (get_shc(l, m).real * get_shc(l2, m-n).real 
                        + get_shc(l, m).imag * get_shc(l2, m-n).imag);
        *l1_sum_imag += get_cg(l1, l2, l, m, n) *
                        (get_shc(l, m).imag * get_shc(l2, m-n).real
                        - get_shc(l, m).real * get_shc(l2, m-n).imag);
    }
}


void calc_l2_sum(const long l1, const long l2, const long l, const long n, double *l2_sum_real, double *l2_sum_imag)
{
    for (long m = maximum(-l, n-l1); m<=minimum(l, n+l1); ++m)
    {
        *l2_sum_real += get_cg(l1, l2, l, m, m-n) *
                        (get_shc(l, m).real * get_shc(l1, m-n).real 
                        + get_shc(l, m).imag * get_shc(l1, m-n).imag);
        *l2_sum_imag += get_cg(l1, l2, l, m, m-n) *
                        (get_shc(l, m).imag * get_shc(l1, m-n).real
                        - get_shc(l, m).real * get_shc(l1, m-n).imag);
    }
}

void bisp_grad()
{
    nnz_ind = 0;

    long col_ind;

    double l1_sum_real, l1_sum_imag;
    double l2_sum_real, l2_sum_imag;
    double l_sum_real, l_sum_imag;

    for (long l1 = 0; l1<=bandlimit; ++l1)
    {
        for (long l2 = 0; l2<=l1; ++l2)
        {
            for (long l = l1-l2; l<=l1+l2 && l<=bandlimit; ++l)
            {
                if (l1==l2 && l==l1)  // if l1=l2=l
                {
                    // Calculate derivative of b_{l1,l2,l} w.r.t real(f_{l,n}) and imag(f_{l,n})
                    for (long n = -l; n<=l; ++n)
                    {
                        col_ind = 2*(l*(l+1) + n);
                        
                        l1_sum_real = 0; 
                        l1_sum_imag = 0;
                        l2_sum_real = 0; 
                        l2_sum_imag = 0;
                        l_sum_real = 0;
                        l_sum_imag = 0;

                        calc_l1_sum(l1, l2, l, n, &l1_sum_real, &l1_sum_imag);
                        calc_l2_sum(l1, l2, l, n, &l2_sum_real, &l2_sum_imag);
                        calc_l_sum(l1, l2, l, n, &l_sum_real, &l_sum_imag);

                         // Save the derivative of real(b_{l1,l2,l}) w.r.t real(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
                        grad_j[nnz_ind] = col_ind;
                        grad_vals[nnz_ind++] = l_sum_real + l1_sum_real + l2_sum_real;

                        // Save the derivative of imag(b_{l1,l2,l}) w.r.t real(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
                        grad_j[nnz_ind] = col_ind;
                        grad_vals[nnz_ind++] = l_sum_imag + l1_sum_imag + l2_sum_imag;

                        // Save the derivative of real(b_{l1,l2,l}) w.r.t imag(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
                        grad_j[nnz_ind] = ++col_ind;
                        grad_vals[nnz_ind++] = -(l_sum_imag - l1_sum_imag - l2_sum_imag);

                        // Save the derivative of imag(b_{l1,l2,l}) w.r.t imag(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
                        grad_j[nnz_ind] = col_ind;
                        grad_vals[nnz_ind++] = l_sum_real - l1_sum_real - l2_sum_real;
                        
                    }
                }
                else if (l1==l2)  // if l1=l2 and l!=l1
                {
                    // Calculate derivative of b_{l1,l2,l} w.r.t real(f_{l,n}) and imag(f_{l,n})
                    for (long n = -l; n<=l; ++n)
                    {
                        col_ind = 2*(l*(l+1) + n);
                        
                        l_sum_real = 0;
                        l_sum_imag = 0;

                        calc_l_sum(l1, l2, l, n, &l_sum_real, &l_sum_imag);

                         // Save the derivative of real(b_{l1,l2,l}) w.r.t real(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
                        grad_j[nnz_ind] = col_ind;
                        grad_vals[nnz_ind++] = l_sum_real;

                        // Save the derivative of imag(b_{l1,l2,l}) w.r.t real(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
                        grad_j[nnz_ind] = col_ind;
                        grad_vals[nnz_ind++] = l_sum_imag;

                        // Save the derivative of real(b_{l1,l2,l}) w.r.t imag(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
                        grad_j[nnz_ind] = ++col_ind;
                        grad_vals[nnz_ind++] = -l_sum_imag;

                        // Save the derivative of imag(b_{l1,l2,l}) w.r.t imag(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
                        grad_j[nnz_ind] = col_ind;
                        grad_vals[nnz_ind++] = l_sum_real;
                        
                    }

                    // Calculate derivative of b_{l1,l2,l} w.r.t real(f_{l1,n}) and imag(f_{l1,n})
                    for (long n = -l1; n<=l1; ++n)
                    {
                        col_ind = 2*(l1*(l1+1) + n);
                        
                        l1_sum_real = 0; 
                        l1_sum_imag = 0;
                        l2_sum_real = 0; 
                        l2_sum_imag = 0;

                        calc_l1_sum(l1, l2, l, n, &l1_sum_real, &l1_sum_imag);
                        calc_l2_sum(l1, l2, l, n, &l2_sum_real, &l2_sum_imag);

                         // Save the derivative of real(b_{l1,l2,l}) w.r.t real(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
                        grad_j[nnz_ind] = col_ind;
                        grad_vals[nnz_ind++] = l1_sum_real + l2_sum_real;

                        // Save the derivative of imag(b_{l1,l2,l}) w.r.t real(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
                        grad_j[nnz_ind] = col_ind;
                        grad_vals[nnz_ind++] = l1_sum_imag + l2_sum_imag;

                        // Save the derivative of real(b_{l1,l2,l}) w.r.t imag(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
                        grad_j[nnz_ind] = ++col_ind;
                        grad_vals[nnz_ind++] = l1_sum_imag + l2_sum_imag;

                        // Save the derivative of imag(b_{l1,l2,l}) w.r.t imag(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
                        grad_j[nnz_ind] = col_ind;
                        grad_vals[nnz_ind++] = - l1_sum_real - l2_sum_real;
                        
                    }
                }
                else if (l==l1)  // if l1=l and l!=l2
                {
                    // Calculate derivative of b_{l1,l2,l} w.r.t real(f_{l2,n}) and imag(f_{l2,n})
                    for (long n = -l2; n<=l2; ++n)
                    {
                        col_ind = 2*(l2*(l2+1) + n);

                        l2_sum_real = 0; 
                        l2_sum_imag = 0;

                        calc_l2_sum(l1, l2, l, n, &l2_sum_real, &l2_sum_imag);

                         // Save the derivative of real(b_{l1,l2,l}) w.r.t real(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
                        grad_j[nnz_ind] = col_ind;
                        grad_vals[nnz_ind++] = l2_sum_real;

                        // Save the derivative of imag(b_{l1,l2,l}) w.r.t real(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
                        grad_j[nnz_ind] = col_ind;
                        grad_vals[nnz_ind++] = l2_sum_imag;

                        // Save the derivative of real(b_{l1,l2,l}) w.r.t imag(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
                        grad_j[nnz_ind] = ++col_ind;
                        grad_vals[nnz_ind++] = l2_sum_imag;

                        // Save the derivative of imag(b_{l1,l2,l}) w.r.t imag(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
                        grad_j[nnz_ind] = col_ind;
                        grad_vals[nnz_ind++] = -l2_sum_real;
                        
                    }

                    // Calculate derivative of b_{l1,l2,l} w.r.t real(f_{l,n}) and imag(f_{l,n})
                    for (long n = -l; n<=l; ++n)
                    {
                        col_ind = 2*(l*(l+1) + n);
                        
                        l1_sum_real = 0; 
                        l1_sum_imag = 0;
                        l_sum_real = 0;
                        l_sum_imag = 0;

                        calc_l1_sum(l1, l2, l, n, &l1_sum_real, &l1_sum_imag);
                        calc_l_sum(l1, l2, l, n, &l_sum_real, &l_sum_imag);

                         // Save the derivative of real(b_{l1,l2,l}) w.r.t real(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
                        grad_j[nnz_ind] = col_ind;
                        grad_vals[nnz_ind++] = l_sum_real + l1_sum_real;

                        // Save the derivative of imag(b_{l1,l2,l}) w.r.t real(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
                        grad_j[nnz_ind] = col_ind;
                        grad_vals[nnz_ind++] = l_sum_imag + l1_sum_imag;

                        // Save the derivative of real(b_{l1,l2,l}) w.r.t imag(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
                        grad_j[nnz_ind] = ++col_ind;
                        grad_vals[nnz_ind++] = -(l_sum_imag - l1_sum_imag);

                        // Save the derivative of imag(b_{l1,l2,l}) w.r.t imag(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
                        grad_j[nnz_ind] = col_ind;
                        grad_vals[nnz_ind++] = l_sum_real - l1_sum_real;
                        
                    }
                }
                else if(l==l2)  // if l2=l and l!=l1
                {
                    // Calculate derivative of b_{l1,l2,l} w.r.t real(f_{l1,n}) and imag(f_{l1,n})
                    for (long n = -l1; n<=l1; ++n)
                    {
                        col_ind = 2*(l1*(l1+1) + n);
                        
                        l1_sum_real = 0; 
                        l1_sum_imag = 0;

                        calc_l1_sum(l1, l2, l, n, &l1_sum_real, &l1_sum_imag);

                         // Save the derivative of real(b_{l1,l2,l}) w.r.t real(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
                        grad_j[nnz_ind] = col_ind;
                        grad_vals[nnz_ind++] = l1_sum_real;

                        // Save the derivative of imag(b_{l1,l2,l}) w.r.t real(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
                        grad_j[nnz_ind] = col_ind;
                        grad_vals[nnz_ind++] = l1_sum_imag;

                        // Save the derivative of real(b_{l1,l2,l}) w.r.t imag(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
                        grad_j[nnz_ind] = ++col_ind;
                        grad_vals[nnz_ind++] = l1_sum_imag;

                        // Save the derivative of imag(b_{l1,l2,l}) w.r.t imag(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
                        grad_j[nnz_ind] = col_ind;
                        grad_vals[nnz_ind++] = - l1_sum_real;
                        
                    }

                    // Calculate derivative of b_{l1,l2,l} w.r.t real(f_{l,n}) and imag(f_{l,n})
                    for (long n = -l; n<=l; ++n)
                    {
                        col_ind = 2*(l*(l+1) + n);
                        
                        l2_sum_real = 0; 
                        l2_sum_imag = 0;
                        l_sum_real = 0;
                        l_sum_imag = 0;

                        calc_l2_sum(l1, l2, l, n, &l2_sum_real, &l2_sum_imag);
                        calc_l_sum(l1, l2, l, n, &l_sum_real, &l_sum_imag);

                         // Save the derivative of real(b_{l1,l2,l}) w.r.t real(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
                        grad_j[nnz_ind] = col_ind;
                        grad_vals[nnz_ind++] = l_sum_real + l2_sum_real;

                        // Save the derivative of imag(b_{l1,l2,l}) w.r.t real(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
                        grad_j[nnz_ind] = col_ind;
                        grad_vals[nnz_ind++] = l_sum_imag + l2_sum_imag;

                        // Save the derivative of real(b_{l1,l2,l}) w.r.t imag(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
                        grad_j[nnz_ind] = ++col_ind;
                        grad_vals[nnz_ind++] = -(l_sum_imag - l2_sum_imag);

                        // Save the derivative of imag(b_{l1,l2,l}) w.r.t imag(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
                        grad_j[nnz_ind] = col_ind;
                        grad_vals[nnz_ind++] = l_sum_real - l2_sum_real;
                        
                    }

                }
                else   // if l1!=l2, l!=l1 and l!=l2
                {
                    // Calculate derivative of b_{l1,l2,l} w.r.t real(f_{l1,n}) and imag(f_{l1,n})
                    for (long n = -l1; n<=l1; ++n)
                    {
                        col_ind = 2*(l1*(l1+1) + n);
                        
                        l1_sum_real = 0; 
                        l1_sum_imag = 0;

                        calc_l1_sum(l1, l2, l, n, &l1_sum_real, &l1_sum_imag);

                         // Save the derivative of real(b_{l1,l2,l}) w.r.t real(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
                        grad_j[nnz_ind] = col_ind;
                        grad_vals[nnz_ind++] = l1_sum_real;

                        // Save the derivative of imag(b_{l1,l2,l}) w.r.t real(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
                        grad_j[nnz_ind] = col_ind;
                        grad_vals[nnz_ind++] = l1_sum_imag;

                        // Save the derivative of real(b_{l1,l2,l}) w.r.t imag(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
                        grad_j[nnz_ind] = ++col_ind;
                        grad_vals[nnz_ind++] = l1_sum_imag;

                        // Save the derivative of imag(b_{l1,l2,l}) w.r.t imag(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
                        grad_j[nnz_ind] = col_ind;
                        grad_vals[nnz_ind++] = -l1_sum_real;
                        
                    }

                    // Calculate derivative of b_{l1,l2,l} w.r.t real(f_{l2,n}) and imag(f_{l2,n})
                    for (long n = -l2; n<=l2; ++n)
                    {
                        col_ind = 2*(l2*(l2+1) + n);
                        
                        l2_sum_real = 0; 
                        l2_sum_imag = 0;

                        calc_l2_sum(l1, l2, l, n, &l2_sum_real, &l2_sum_imag);

                         // Save the derivative of real(b_{l1,l2,l}) w.r.t real(f_{l2,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
                        grad_j[nnz_ind] = col_ind;
                        grad_vals[nnz_ind++] =  l2_sum_real;

                        // Save the derivative of imag(b_{l1,l2,l}) w.r.t real(f_{l2,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
                        grad_j[nnz_ind] = col_ind;
                        grad_vals[nnz_ind++] = l2_sum_imag;

                        // Save the derivative of real(b_{l1,l2,l}) w.r.t imag(f_{l2,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
                        grad_j[nnz_ind] = ++col_ind;
                        grad_vals[nnz_ind++] = l2_sum_imag;

                        // Save the derivative of imag(b_{l1,l2,l}) w.r.t imag(f_{l2,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
                        grad_j[nnz_ind] = col_ind;
                        grad_vals[nnz_ind++] = -l2_sum_real;
                        
                    }

                    // Calculate derivative of b_{l1,l2,l} w.r.t real(f_{l,n}) and imag(f_{l,n})
                    for (long n = -l; n<=l; ++n)
                    {
                        col_ind = 2*(l*(l+1) + n);
                        
                        l_sum_real = 0;
                        l_sum_imag = 0;

                        calc_l_sum(l1, l2, l, n, &l_sum_real, &l_sum_imag);

                         // Save the derivative of real(b_{l1,l2,l}) w.r.t real(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
                        grad_j[nnz_ind] = col_ind;
                        grad_vals[nnz_ind++] = l_sum_real;

                        // Save the derivative of imag(b_{l1,l2,l}) w.r.t real(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
                        grad_j[nnz_ind] = col_ind;
                        grad_vals[nnz_ind++] = l_sum_imag;

                        // Save the derivative of real(b_{l1,l2,l}) w.r.t imag(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
                        grad_j[nnz_ind] = ++col_ind;
                        grad_vals[nnz_ind++] = -l_sum_imag;

                        // Save the derivative of imag(b_{l1,l2,l}) w.r.t imag(f_{l,n})
                        grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
                        grad_j[nnz_ind] = col_ind;
                        grad_vals[nnz_ind++] = l_sum_real;
                        
                    }

                }

            }
        }
    }
}

// void bisp_grad()
// {
//     nnz_ind = 0;
    
//     long row_ind;
//     long col_ind;

//     long l1;
//     long l2;
//     long l;

//     double der_real;
//     double der_imag;
//     double der_real2;
//     double der_imag2;
//     double der_real3;
//     double der_imag3;

//     // Case 1: Calculate the derivative of b_{l1,l2,l} for l1<l2, l!=l1 and l!=l2
//     for (l1 = 0; l1<=bandlimit; ++l1)
//     {
//         for (l2 = 0; l2<l1; ++l2)
//         {
//             // Calculate the derivative of b_{l1,l2,l} for l1<l2, l<l2
//             for (l = l1-l2; l<l2; ++l)
//             {
//                 // Derivative w.r.t real(f_{l,m}) and imag(f_{l,m})
//                 for (long m = -l; m<=l; ++m)
//                 {
//                     der_real = 0;
//                     der_imag = 0;

//                     col_ind = 2*(l*(l+1) + m); // Column index of the derivative w.r.t to real(f_{l,m})

//                     for (long m1 = maximum(-l1, m-l2); m1<=minimum(l1, m+l2); ++m1)
//                     {
//                         der_real += get_cg(l1, l2, l, m, m1) 
//                                     * (get_shc(l1, m1).real * get_shc(l2, m-m1).real 
//                                         - get_shc(l1, m1).imag * get_shc(l2, m-m1).imag);
//                         der_imag += - get_cg(l1, l2, l, m, m1) 
//                                     * (get_shc(l1, m1).imag * get_shc(l2, m-m1).real 
//                                         + get_shc(l1, m1).real * get_shc(l2, m-m1).imag);
//                     }

//                     // Save the derivative of real(b_{l1,l2,l}) w.r.t real(f_{l,m})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
//                     grad_j[nnz_ind] = col_ind;
//                     grad_vals[nnz_ind++] = der_real;

//                     // Save the derivative of imag(b_{l1,l2,l}) w.r.t real(f_{l,m})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
//                     grad_j[nnz_ind] = col_ind;
//                     grad_vals[nnz_ind++] = der_imag;

//                     // Save the derivative of real(b_{l1,l2,l}) w.r.t imag(f_{l,m})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
//                     grad_j[nnz_ind] = ++col_ind;
//                     grad_vals[nnz_ind++] = -der_imag;

//                     // Save the derivative of imag(b_{l1,l2,l}) w.r.t imag(f_{l,m})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
//                     grad_j[nnz_ind] = col_ind;
//                     grad_vals[nnz_ind++] = der_real;
                    
//                 }


//                 // Derivative w.r.t real(f_{l1,m1}) and imag(f_{l1,m1})
//                 for (long m1 = -l1; m1<=l1; ++m1)
//                 {
//                     der_real = 0;
//                     der_imag = 0;

//                     col_ind = 2*(l1*(l1+1) + m1); // Column index of the derivative w.r.t to real(f_{l1,m1})

//                     for (long m = maximum(-l, m1-l2); m<=minimum(l, m1+l2); ++m)
//                     {
//                         der_real += get_cg(l1, l2, l, m, m-m1) 
//                                     * (get_shc(l, m).real * get_shc(l2, m-m1).real 
//                                         + get_shc(l, m).imag * get_shc(l2, m-m1).imag);
//                         der_imag += get_cg(l1, l2, l, m, m-m1) 
//                                     * (get_shc(l, m).imag * get_shc(l2, m-m1).real 
//                                         + get_shc(l, m).real * get_shc(l2, m-m1).imag);
//                     }

//                     // Save the derivative of real(b_{l1,l2,l}) w.r.t real(f_{l1,m1})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
//                     grad_j[nnz_ind] = col_ind;
//                     grad_vals[nnz_ind++] = der_real;

//                     // Save the derivative of imag(b_{l1,l2,l}) w.r.t real(f_{l1,m1})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
//                     grad_j[nnz_ind] = col_ind;
//                     grad_vals[nnz_ind++] = der_imag;

//                     // Save the derivative of real(b_{l1,l2,l}) w.r.t imag(f_{l1,m1})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
//                     grad_j[nnz_ind] = ++col_ind;
//                     grad_vals[nnz_ind++] = -der_imag;

//                     // Save the derivative of imag(b_{l1,l2,l}) w.r.t imag(f_{l1,m1})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
//                     grad_j[nnz_ind] = col_ind;
//                     grad_vals[nnz_ind++] = -der_real;
//                 }


//                 // Derivative w.r.t real(f_{l2,m2}) and imag(f_{l2,m2})
//                 for (long m2 = -l2; m2<=l2; ++m2)
//                 {
//                     der_real = 0;
//                     der_imag = 0;

//                     col_ind = 2*(l2*(l2+1) + m2); // Column index of the derivative w.r.t to real(f_{l2,m2})

//                     for (long m = maximum(-l, m2-l2); m<=minimum(l, m2+l2); ++m)
//                     {
//                         der_real += get_cg(l1, l2, l, m, m-m2) 
//                                     * (get_shc(l, m).real * get_shc(l1, m-m2).real 
//                                         + get_shc(l, m).imag * get_shc(l1, m-m2).imag);
//                         der_imag += get_cg(l1, l2, l, m, m-m2) 
//                                     * (get_shc(l, m).imag * get_shc(l1, m-m2).real 
//                                         + get_shc(l, m).real * get_shc(l1, m-m2).imag);
//                     }

//                     // Save the derivative of real(b_{l1,l2,l}) w.r.t real(f_{l2,m2})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
//                     grad_j[nnz_ind] = col_ind;
//                     grad_vals[nnz_ind++] = der_real;

//                     // Save the derivative of imag(b_{l1,l2,l}) w.r.t real(f_{l2,m2})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
//                     grad_j[nnz_ind] = col_ind;
//                     grad_vals[nnz_ind++] = der_imag;

//                     // Save the derivative of real(b_{l1,l2,l}) w.r.t imag(f_{l2,m2})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
//                     grad_j[nnz_ind] = ++col_ind;
//                     grad_vals[nnz_ind++] = -der_imag;

//                     // Save the derivative of imag(b_{l1,l2,l}) w.r.t imag(f_{l2,m2})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
//                     grad_j[nnz_ind] = col_ind;
//                     grad_vals[nnz_ind++] = -der_real;
//                 }
//             }

//             // Calculate the derivative of b_{l1,l2,l} for l1<l2, l2<l<l1
//             for (l = l2+1; l<l1; ++l)
//             {
//                 // Derivative w.r.t real(f_{l,m}) and imag(f_{l,m})
//                 for (long m = -l; m<=l; ++m)
//                 {
//                     der_real = 0;
//                     der_imag = 0;

//                     col_ind = 2*(l*(l+1) + m); // Column index of the derivative w.r.t to real(f_{l,m})

//                     for (long m1 = maximum(-l1, m-l2); m1<=minimum(l1, m+l2); ++m1)
//                     {
//                         der_real += get_cg(l1, l2, l, m, m1) 
//                                     * (get_shc(l1, m1).real * get_shc(l2, m-m1).real 
//                                         - get_shc(l1, m1).imag * get_shc(l2, m-m1).imag);
//                         der_imag += - get_cg(l1, l2, l, m, m1) 
//                                     * (get_shc(l1, m1).imag * get_shc(l2, m-m1).real 
//                                         + get_shc(l1, m1).real * get_shc(l2, m-m1).imag);
//                     }

//                     // Save the derivative of real(b_{l1,l2,l}) w.r.t real(f_{l,m})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
//                     grad_j[nnz_ind] = col_ind;
//                     grad_vals[nnz_ind++] = der_real;

//                     // Save the derivative of imag(b_{l1,l2,l}) w.r.t real(f_{l,m})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
//                     grad_j[nnz_ind] = col_ind;
//                     grad_vals[nnz_ind++] = der_imag;

//                     // Save the derivative of real(b_{l1,l2,l}) w.r.t imag(f_{l,m})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
//                     grad_j[nnz_ind] = ++col_ind;
//                     grad_vals[nnz_ind++] = -der_imag;

//                     // Save the derivative of imag(b_{l1,l2,l}) w.r.t imag(f_{l,m})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
//                     grad_j[nnz_ind] = col_ind;
//                     grad_vals[nnz_ind++] = der_real;
                    
//                 }


//                 // Derivative w.r.t real(f_{l1,m1}) and imag(f_{l1,m1})
//                 for (long m1 = -l1; m1<=l1; ++m1)
//                 {
//                     der_real = 0;
//                     der_imag = 0;

//                     col_ind = 2*(l1*(l1+1) + m1); // Column index of the derivative w.r.t to real(f_{l1,m1})

//                     for (long m = maximum(-l, m1-l2); m<=minimum(l, m1+l2); ++m)
//                     {
//                         der_real += get_cg(l1, l2, l, m, m-m1) 
//                                     * (get_shc(l, m).real * get_shc(l2, m-m1).real 
//                                         + get_shc(l, m).imag * get_shc(l2, m-m1).imag);
//                         der_imag += get_cg(l1, l2, l, m, m-m1) 
//                                     * (get_shc(l, m).imag * get_shc(l2, m-m1).real 
//                                         + get_shc(l, m).real * get_shc(l2, m-m1).imag);
//                     }

//                     // Save the derivative of real(b_{l1,l2,l}) w.r.t real(f_{l1,m1})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
//                     grad_j[nnz_ind] = col_ind;
//                     grad_vals[nnz_ind++] = der_real;

//                     // Save the derivative of imag(b_{l1,l2,l}) w.r.t real(f_{l1,m1})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
//                     grad_j[nnz_ind] = col_ind;
//                     grad_vals[nnz_ind++] = der_imag;

//                     // Save the derivative of real(b_{l1,l2,l}) w.r.t imag(f_{l1,m1})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
//                     grad_j[nnz_ind] = ++col_ind;
//                     grad_vals[nnz_ind++] = -der_imag;

//                     // Save the derivative of imag(b_{l1,l2,l}) w.r.t imag(f_{l1,m1})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
//                     grad_j[nnz_ind] = col_ind;
//                     grad_vals[nnz_ind++] = -der_real;
//                 }


//                 // Derivative w.r.t real(f_{l2,m2}) and imag(f_{l2,m2})
//                 for (long m2 = -l2; m2<=l2; ++m2)
//                 {
//                     der_real = 0;
//                     der_imag = 0;

//                     col_ind = 2*(l2*(l2+1) + m2); // Column index of the derivative w.r.t to real(f_{l2,m2})

//                     for (long m = maximum(-l, m2-l2); m<=minimum(l, m2+l2); ++m)
//                     {
//                         der_real += get_cg(l1, l2, l, m, m-m2) 
//                                     * (get_shc(l, m).real * get_shc(l1, m-m2).real 
//                                         + get_shc(l, m).imag * get_shc(l1, m-m2).imag);
//                         der_imag += get_cg(l1, l2, l, m, m-m2) 
//                                     * (get_shc(l, m).imag * get_shc(l1, m-m2).real 
//                                         + get_shc(l, m).real * get_shc(l1, m-m2).imag);
//                     }

//                     // Save the derivative of real(b_{l1,l2,l}) w.r.t real(f_{l2,m2})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
//                     grad_j[nnz_ind] = col_ind;
//                     grad_vals[nnz_ind++] = der_real;

//                     // Save the derivative of imag(b_{l1,l2,l}) w.r.t real(f_{l2,m2})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
//                     grad_j[nnz_ind] = col_ind;
//                     grad_vals[nnz_ind++] = der_imag;

//                     // Save the derivative of real(b_{l1,l2,l}) w.r.t imag(f_{l2,m2})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
//                     grad_j[nnz_ind] = ++col_ind;
//                     grad_vals[nnz_ind++] = -der_imag;

//                     // Save the derivative of imag(b_{l1,l2,l}) w.r.t imag(f_{l2,m2})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
//                     grad_j[nnz_ind] = col_ind;
//                     grad_vals[nnz_ind++] = -der_real;
//                 }
//             }

//             // Calculate the derivative of b_{l1,l2,l} for l1<l2, l<=l1+l2
//             for (l = l1+1; l<=l1+l2 && l<=bandlimit; ++l)
//             {
//                 // Derivative w.r.t real(f_{l,m}) and imag(f_{l,m})
//                 for (long m = -l; m<=l; ++m)
//                 {
//                     der_real = 0;
//                     der_imag = 0;

//                     col_ind = 2*(l*(l+1) + m); // Column index of the derivative w.r.t to real(f_{l,m})

//                     for (long m1 = maximum(-l1, m-l2); m1<=minimum(l1, m+l2); ++m1)
//                     {
//                         der_real += get_cg(l1, l2, l, m, m1) 
//                                     * (get_shc(l1, m1).real * get_shc(l2, m-m1).real 
//                                         - get_shc(l1, m1).imag * get_shc(l2, m-m1).imag);
//                         der_imag += - get_cg(l1, l2, l, m, m1) 
//                                     * (get_shc(l1, m1).imag * get_shc(l2, m-m1).real 
//                                         + get_shc(l1, m1).real * get_shc(l2, m-m1).imag);
//                     }

//                     // Save the derivative of real(b_{l1,l2,l}) w.r.t real(f_{l,m})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
//                     grad_j[nnz_ind] = col_ind;
//                     grad_vals[nnz_ind++] = der_real;

//                     // Save the derivative of imag(b_{l1,l2,l}) w.r.t real(f_{l,m})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
//                     grad_j[nnz_ind] = col_ind;
//                     grad_vals[nnz_ind++] = der_imag;

//                     // Save the derivative of real(b_{l1,l2,l}) w.r.t imag(f_{l,m})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
//                     grad_j[nnz_ind] = ++col_ind;
//                     grad_vals[nnz_ind++] = -der_imag;

//                     // Save the derivative of imag(b_{l1,l2,l}) w.r.t imag(f_{l,m})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
//                     grad_j[nnz_ind] = col_ind;
//                     grad_vals[nnz_ind++] = der_real;
                    
//                 }


//                 // Derivative w.r.t real(f_{l1,m1}) and imag(f_{l1,m1})
//                 for (long m1 = -l1; m1<=l1; ++m1)
//                 {
//                     der_real = 0;
//                     der_imag = 0;

//                     col_ind = 2*(l1*(l1+1) + m1); // Column index of the derivative w.r.t to real(f_{l1,m1})

//                     for (long m = maximum(-l, m1-l2); m<=minimum(l, m1+l2); ++m)
//                     {
//                         der_real += get_cg(l1, l2, l, m, m-m1) 
//                                     * (get_shc(l, m).real * get_shc(l2, m-m1).real 
//                                         + get_shc(l, m).imag * get_shc(l2, m-m1).imag);
//                         der_imag += get_cg(l1, l2, l, m, m-m1) 
//                                     * (get_shc(l, m).imag * get_shc(l2, m-m1).real 
//                                         + get_shc(l, m).real * get_shc(l2, m-m1).imag);
//                     }

//                     // Save the derivative of real(b_{l1,l2,l}) w.r.t real(f_{l1,m1})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
//                     grad_j[nnz_ind] = col_ind;
//                     grad_vals[nnz_ind++] = der_real;

//                     // Save the derivative of imag(b_{l1,l2,l}) w.r.t real(f_{l1,m1})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
//                     grad_j[nnz_ind] = col_ind;
//                     grad_vals[nnz_ind++] = der_imag;

//                     // Save the derivative of real(b_{l1,l2,l}) w.r.t imag(f_{l1,m1})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
//                     grad_j[nnz_ind] = ++col_ind;
//                     grad_vals[nnz_ind++] = -der_imag;

//                     // Save the derivative of imag(b_{l1,l2,l}) w.r.t imag(f_{l1,m1})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
//                     grad_j[nnz_ind] = col_ind;
//                     grad_vals[nnz_ind++] = -der_real;
//                 }


//                 // Derivative w.r.t real(f_{l2,m2}) and imag(f_{l2,m2})
//                 for (long m2 = -l2; m2<=l2; ++m2)
//                 {
//                     der_real = 0;
//                     der_imag = 0;

//                     col_ind = 2*(l2*(l2+1) + m2); // Column index of the derivative w.r.t to real(f_{l2,m2})

//                     for (long m = maximum(-l, m2-l2); m<=minimum(l, m2+l2); ++m)
//                     {
//                         der_real += get_cg(l1, l2, l, m, m-m2) 
//                                     * (get_shc(l, m).real * get_shc(l1, m-m2).real 
//                                         + get_shc(l, m).imag * get_shc(l1, m-m2).imag);
//                         der_imag += get_cg(l1, l2, l, m, m-m2) 
//                                     * (get_shc(l, m).imag * get_shc(l1, m-m2).real 
//                                         + get_shc(l, m).real * get_shc(l1, m-m2).imag);
//                     }

//                     // Save the derivative of real(b_{l1,l2,l}) w.r.t real(f_{l2,m2})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
//                     grad_j[nnz_ind] = col_ind;
//                     grad_vals[nnz_ind++] = der_real;

//                     // Save the derivative of imag(b_{l1,l2,l}) w.r.t real(f_{l2,m2})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
//                     grad_j[nnz_ind] = col_ind;
//                     grad_vals[nnz_ind++] = der_imag;

//                     // Save the derivative of real(b_{l1,l2,l}) w.r.t imag(f_{l2,m2})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
//                     grad_j[nnz_ind] = ++col_ind;
//                     grad_vals[nnz_ind++] = -der_imag;

//                     // Save the derivative of imag(b_{l1,l2,l}) w.r.t imag(f_{l2,m2})
//                     grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
//                     grad_j[nnz_ind] = col_ind;
//                     grad_vals[nnz_ind++] = -der_real;
//                 }
//             }
//         }
//     }

//     // Case 2: Calculate the derivative of b_{l1,l2,l} for l1<l2, l=l1 or l=l2
//     for (l1 = 0; l1<=bandlimit; ++l1)
//     {
//         for (l2 = 0; l2<l1; ++l2)
//         {
//             // Calculate the derivative of b_{l1,l2,l} for l1<l2, l=l1
//             l = l1;

//             // Derivative w.r.t real(f_{l,m}) and imag(f_{l,m})
//             for (long m = -l; m<=l; ++m)
//             {
//                 der_real = 0;
//                 der_imag = 0;

//                 col_ind = 2*(l*(l+1) + m); // Column index of the derivative w.r.t to real(f_{l,m})

//                 for (long m1 = maximum(-l1, m-l2); m1<=minimum(l1, m+l2); ++m1)
//                 {
//                     der_real += get_cg(l1, l2, l, m, m-m1) 
//                                 * (get_shc(l1, m1).real * get_shc(l2, m-m1).real
//                                     - get_shc(l1, m1).imag * get_shc(l2, m-m1).imag);
//                     der_imag += - get_cg(l1, l2, l, m, m-m1) *
//                                 (get_shc(l1, m1).imag * get_shc(l2, m-m1).real 
//                                 + get_shc(l1, m1).real * get_shc(l2, m-m1).imag) ;
//                 }

//                 // Save the derivative of real(b_{l1,l2,l}) w.r.t real(f_{l,m})
//                 grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
//                 grad_j[nnz_ind] = col_ind;
//                 grad_vals[nnz_ind++] = der_real;

//                 // Save the derivative of imag(b_{l1,l2,l}) w.r.t real(f_{l,m})
//                 grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
//                 grad_j[nnz_ind] = col_ind;
//                 grad_vals[nnz_ind++] = der_imag;

//                 // Save the derivative of real(b_{l1,l2,l}) w.r.t imag(f_{l,m})
//                 grad_i[nnz_ind] = bisp_lookup(l1, l2, l, REAL_PART);
//                 grad_j[nnz_ind] = ++col_ind;
//                 grad_vals[nnz_ind++] = -der_imag;

//                 // Save the derivative of imag(b_{l1,l2,l}) w.r.t imag(f_{l,m})
//                 grad_i[nnz_ind] = bisp_lookup(l1, l2, l, IMAG_PART);
//                 grad_j[nnz_ind] = col_ind;
//                 grad_vals[nnz_ind++] = der_real;
                
//             }
//             // TODO

//             // Derivative w.r.t real(f_{l2,m2}) and imag(f_{l2,m2})
//             // TODO
            

//             // Calculate the derivative of b_{l1,l2,l} for l1<l2, l=l2
//             if (2*l2>=l1)
//             {
//                 l = l2;
//                 // Derivative w.r.t real(f_{l,m}) and imag(f_{l,m})
//                 // TODO

//                 // Derivative w.r.t real(f_{l1,m1}) and imag(f_{l1,m1})
//                 // TODO
//             }
//         }
//     }

//     // Case 3: Calculate the derivative of b_{l1,l2,l} for l1=l2, l!=l1
//     for (l1 = 0; l1<=bandlimit; ++l1)
//     {
//         l2 = l1;

//         // Calculate the derivative of b_{l1,l2,l} for l1=l2, l<l1
//         for (l = 0; l<l1; ++l)
//         {
//             // TODO

//             ++nnz_ind;
//         }

//         // Calculate the derivative of b_{l1,l2,l} for l1=l2, l1<l<l1+l2
//         for (l = l1+1; l<=2*l1 && l<=bandlimit; ++l)
//         {
//             // TODO

//             ++nnz_ind;
//         }
//     }

//     // Case 4: Calculate the derivative of b_{l1,l2,l} for l1=l2=l
//     for (l1 = 0; l1<=bandlimit; ++l1)
//     {
//         l2 = l1;
//         l = l1;

//         // TODO

//         ++nnz_ind;
//     }
// }


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
        mxArray *inputMxArray[1];
        inputMxArray[0] = mxCreateDoubleScalar(bandlimit);
        mexCallMATLAB(0, NULL, 1, inputMxArray, "loadCGTable");

        CGs = mexGetVariablePtr("global", "CGs");

        destory_bisp_lookup_table();
        destroy_CGTable();

        previous_bandlimit = bandlimit;

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
        bisp_grad();
    }
}
