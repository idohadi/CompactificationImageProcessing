/** 
 * Apply a rotation in unit quaternion representation to points.
 * 
 * MATLAB call format:
 *      rotatedPoints = applyRotation_mex(quaternion, points)
 * where
 *  quaternion is a 4 x 1 or 1 x 4 unit array.
 *  points is an 3 x N array of real doubles.
 *  rotatedPoints are 3 x N array such that 
 *    rotatedPoints(:, J) = points(:, J) rotated using the rotation represntated by quaternion.
 * 
 * NOTE:
 *  Code performs no input checks.
  */

#include "mex.h"


double *quaternion;
double *points;
size_t N;

double *rotatedPoints;


void cross_product(double *a, double *b, double *result)
{
    /* Calculate the cross product a x b and save the reuslt in result.
        Code enables result to be identical to a or b. */
    double i = a[1]*b[2] - a[2]*b[1];
    double j = a[2]*b[0] - a[0]*b[2];
    double k = a[0]*b[1] - a[1]*b[0];
    result[0] = i;
    result[1] = j;
    result[2] = k;
}

void apply_rotation(double *x, double *quaternion, double *result)
{
    /* Apply the rotation represented by quternion to x and save the result in result.  */
    cross_product(&quaternion[1], x, result);
    for (int i = 0; i < 3; ++i)
    {
        result[i] += quaternion[0]*x[i];
    }
    cross_product(&quaternion[1], result, result);
    for (int i = 0; i < 3; ++i)
    {
        result[i] *= 2;
        result[i] += x[i];
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Get input
    quaternion = mxGetDoubles(prhs[0]);
    points = mxGetDoubles(prhs[1]);
    N = mxGetN(prhs[1]);

    // Create output
    plhs[0] = mxCreateDoubleMatrix(3, N, mxREAL);
    
    // Calculate output
    for (size_t n = 0; n<N; ++n)
    {
      apply_rotation(&points[3*n], quaternion, &rotatedPoints[3*n]);
    }
}
