// TODO: write docs in doxygen
/** 
 * Evaluate real-valued spherical function, represented by spherical harmonics coefficients.
 * 
 * MATLAB call form:
 *      Yt = evalYt_mex(theta, phi, t)
 *  where
 *      theta and phi and arrays of size N x 1 or 1 x N representing spherical points
 *      t is a non-negative integer representing the maximum order to calculate
 *      Yt is an N x ((t+1)*(t+1)) matrix such that
 *          Yt(n, M) = Y_{l}^{m} (theta(n),phi(n)),
 *      where M is the index of (l,m) in the lexicographical ordering of { (l,m) in Z^2 | 0<=l<=t, abs(m)<=l }
 *      and Y_{l}^{m} is the spherical harmonics of order l and degree m.
 * 
 * NOTE:
 *  (1) Code performs no input checks.
 *  (2) This function was originally written using pure MATLAB in my old codebase.
 */

#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>
#include "mex.h"
#include "spherical_harmonics.h"
#include "alegendre.h"

double *theta;
double *phi;
size_t bandlimit;
size_t N;
bool first_run = true;

mxComplexDouble *Yt;
size_t rows_no;
long ind;

double cos_recur;
double sin_recur;
double cos_val;
double sin_val;
double cos_recur_next;
double sin_recur_next;


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (first_run)
    {
        alegendre_init();
        first_run = false;
        mexPrintf("Algendere succesfully initialized.\n");
    }

    // Get input variables
    theta = mxGetDoubles(prhs[0]);
    phi = mxGetDoubles(prhs[1]);
    bandlimit = (size_t) mxGetScalar(prhs[2]);
    
    N = mxGetNumberOfElements(prhs[0]);

    // Create output
    rows_no = (bandlimit+1)*(bandlimit+1);
    plhs[0] = mxCreateDoubleMatrix(rows_no, N, mxCOMPLEX);
    Yt = mxGetComplexDoubles(plhs[0]);

    #pragma omp parallel for num_threads(10)
    for (long n = 0; n<N; ++n)
    {
        for (long l = 0; l<=bandlimit; ++l)
        {
            cos_recur = cos(phi[n]);
            cos_val = cos_recur;

            sin_recur = sin(phi[n]);
            sin_val = sin_recur;
            
            ind = rows_no*n + l*l+l;

            // Calculate the real part and imaginary part of Y_{l}^{0} (theta(n),phi(n))
            Yt[ind].real = alegendre(l, 0, theta[n]);

            for (long m = 1; m<=l; ++m)
            {
                // Calculate the real part and imaginary part of Y_{l}^{m} (theta(n),phi(n)) and Y_{l}^{-m} (theta(n),phi(n)) 
                Yt[ind+m].real = ((m%2)==0 ? 1 : -1)*alegendre(l, m, theta[n]);
                Yt[ind+m].imag = sin_recur*Yt[ind+m].real;

                Yt[ind-m].real = ((m%2)==0 ? 1 : -1) * Yt[ind+m].real;
                Yt[ind-m].imag = -sin_recur*Yt[ind-m].real;

                Yt[ind-m].real *= cos_recur;
                Yt[ind+m].real *= cos_recur;


                // Calculate cos((n+1)theta[n]) and sin((n+1)theta[n])
                cos_recur_next = cos_val*cos_recur - sin_val*sin_recur;
                sin_recur_next = sin_val*cos_recur + cos_val*sin_recur;

                cos_recur = cos_recur_next;
                sin_recur = sin_recur_next;
            }
        }
    }
}