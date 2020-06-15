// TODO: write docs in doxygen
/** 
 * Return the spherical harmonics coefficients of the composition of a bandlimited spherical function and a rotation (in quaternion representation)
 * 
 * MATLAB call form:
 *      rotated_shc = rotateSHC(shc, bandlimit, tdesign_bandlimit, rotation, real_vs_imag)
 *  where
 *      shc is a column array of size (bandlimit+1)^2 if real_vs_imag=0 and of size 2*(bandlimit+1)^2 if real_vs_imag~=0
 *      rotation is the quaternion representation of a rotation.
 *      tdesign_bandlimit is the bandlimit of the t-design to be used.
 *      real_vs_imag is 0, if shc represents real-valued spherical function and non-zero if it represents a complex-valued spherical function.
 * 
 * 
 * NOTE:
 *  Code performs no input checks.
  */

#include <stdint.h>
#include <stdbool.h>
#include "mex.h"
#include "spherical_harmonics.h"

size_t bandlimit;
size_t previous_tdesign_bandlimit, tdesign_bandlimit;
size_t real_vs_imag;
r_shc shc;

double *rotation;

tdesign_cart td;

bool first_run = true;

r_shc output_shc;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Get input variables
    bandlimit = (size_t) mxGetScalar(prhs[1]);
    tdesign_bandlimit = (size_t) mxGetScalar(prhs[2]);
    rotation = mxGetDoubles(prhs[3]); 
    real_vs_imag = (size_t) mxGetScalar(prhs[4]);

    // Initilize SFMT on the first run
    if (first_run==true)
    {
        td = read_tdesign(tdesign_bandlimit);
        previous_tdesign_bandlimit = tdesign_bandlimit;
        first_run = false;
        mexPrintf("t-design succesfully initialized.\n");
        alegendre_init();
        mexPrintf("Algendere succesfully initialized.\n");
    }

    if (previous_tdesign_bandlimit!=tdesign_bandlimit)
    {
        previous_tdesign_bandlimit = tdesign_bandlimit;
        deallocate_tdesign(&td);
        td = read_tdesign(tdesign_bandlimit);
        mexPrintf("t-design succesfully reinitialized.\n");
    }

    // Create output
    if (real_vs_imag==0)
    {
        plhs[0] = mxCreateDoubleMatrix((bandlimit+1)*(bandlimit+1), 1, mxREAL);
        shc = r_init_shc(bandlimit, mxGetDoubles(prhs[0]));
        output_shc = r_init_shc(bandlimit, mxGetDoubles(plhs[0]));
        r_rotate_spherical_harmonics(&shc, rotation, &td, &output_shc);
    }
    else
    {
        // TODO: handle complex-valued functions in future
    }
}