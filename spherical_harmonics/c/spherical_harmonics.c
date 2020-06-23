// TODO: write docs in doxygen

#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "spherical_harmonics.h"
#include "SFMT.h"
#include "alegendre.h"
#include "utility_functions.h"

#define PI 3.14159265358979323846


double *r_allocate_coefficients(const size_t bandlimit)
{
    /*
    Allocate memory for coefficients.
     */
    
    return malloc((bandlimit+1)*(bandlimit+1)*sizeof(double));
}


double *c_allocate_coefficients(const size_t bandlimit)
{
    /*
    Allocate memory for coefficients.
     */
    
    return malloc(2*(bandlimit+1)*(bandlimit+1)*sizeof(double));
}


void r_deallocate_coefficients(double * const coefficients)
{
    free(coefficients);
}


void c_deallocate_coefficients(double * const coefficients)
{
    free(coefficients);
}


r_shc r_init_shc(const size_t bandlimit, double * const coefficients)
{
    /*
    Create and initialize an r_shc.
    Coefficients are given as pointer to array of appropriate size.
     */
    
    r_shc shc;
    shc.bandlimit = bandlimit;
    shc.coefficients = coefficients;
    return shc;
}


c_shc c_init_shc(const size_t bandlimit, double * const coefficients)
{
    /*
    Create and initialize an c_shc.
    Coefficients are given as pointer to array of appropriate size.
     */
    
    c_shc shc;
    shc.bandlimit = bandlimit;
    shc.coefficients = coefficients;
    return shc;
}


r_shc r_init_shc2(const size_t bandlimit)
{
    /* 
    Initialize an r_shc with the bandlimit and an all-zero array of coefficients. 
    */
    
    r_shc shc;
    shc.bandlimit = bandlimit;
    shc.coefficients = r_allocate_coefficients(bandlimit);
    for (long l = 0; l<=bandlimit; ++l)
    {
        r_set_shc(&shc, REAL_PART, l, 0, 0.0);
        for (long m = 1; m<=l; ++m)
        {
            r_set_shc(&shc, REAL_PART, l, m, 0.0);
            r_set_shc(&shc, IMAG_PART, l, m, 0.0);
        }
    }
    return shc;
}


c_shc c_init_shc2(const size_t bandlimit)
{
    /* 
    Initialize an r_shc with the bandlimit and an all-zero array of coefficients. 
    */
    
    c_shc shc;
    shc.bandlimit = bandlimit;
    shc.coefficients = c_allocate_coefficients(bandlimit);
    for (long l = 0; l<=bandlimit; ++l)
    {
        for (long m = -l; m<=l; ++m)
        {
            c_set_shc(&shc, REAL_PART, l, m, 0.0);
            c_set_shc(&shc, IMAG_PART, l, m, 0.0);
        }
    }
    return shc;
}


double r_get_shc(r_shc const *shc, const PART part, const long l, const long m)
{
   if (l > shc->bandlimit)
   {
       return 0.0;
   }
   else
   {
        if (m==0)
        {
            if (part==IMAG_PART)
            {
                return 0.0;
            }
            else if (part==REAL_PART)
            {
                return shc->coefficients[l*l];
            }
        }

        if (m>0)
        {
            if (part==REAL_PART)
            {
                return shc->coefficients[l*l + 2*m - 1];
            }
            else if (part==IMAG_PART)
            {
                return shc->coefficients[l*l + 2*m];
            }
        }

        if (m<0)
        {
            if (part==REAL_PART)
            {
                return ((m % 2 == 0) ? 1 : -1)*shc->coefficients[l*l - 2*m - 1];
            }
            else if (part==IMAG_PART)
            {
                return ((m % 2 == 0) ? -1 : 1)*shc->coefficients[l*l - 2*m];
            }
        }
    }
}


double c_get_shc(c_shc const *shc, const PART part, const long l, const long m)
{
    if (l > shc->bandlimit)
    {
        return 0.0;
    }
    else
    {
        if (part == REAL_PART)
        {
            return shc->coefficients[2*l*l + 2*(m + l) ];
        }
        else if (part == IMAG_PART)
        {
            return shc->coefficients[2*l*l + 2*(m + l) + 1];
        }
    }
}


bool r_set_shc(r_shc const *shc, const PART part, const long l, const long m, const double value)
{
    if (l > shc->bandlimit)
   {
       return false;
   }
   else
   {
        if (m==0)
        {
            if (part==IMAG_PART)
            {
                return false;
            }
            else if (part==REAL_PART)
            {
                shc->coefficients[l*l] = value;
                return true;
            }
        }

        if (m>0)
        {
            if (part==REAL_PART)
            {
                shc->coefficients[l*l + 2*m - 1] = value;
                return true;
            }
            else if (part==IMAG_PART)
            {
                shc->coefficients[l*l + 2*m] = value;
                return true;
            }
        }

        if (m<0)
        {
            if (part==REAL_PART)
            {
                shc->coefficients[l*l - 2*m - 1] = ((m % 2 == 0) ? 1 : -1)*value;
                return true;
            }
            else if (part==IMAG_PART)
            {
                shc->coefficients[l*l - 2*m] = ((m % 2 == 0) ? -1 : 1)*value;
                return true;
            }
        }
    }
}

long r_lm_to_index(const PART part, const long l, const long m)
{
    if (m==0)
    {
        return l*l;
    }

    if (m>0)
    {
        if (part==REAL_PART)
        {
            return l*l + 2*m - 1;
        }
        else if (part==IMAG_PART)
        {
            return l*l + 2*m;
        }
    }

    if (m<0)
    {
        if (part==REAL_PART)
        {
            return l*l - 2*m - 1;
        }
        else if (part==IMAG_PART)
        {
            return l*l - 2*m;
        }
    }
}

long c_lm_to_index(const PART part, const long l, const long m)
{
    if (part == REAL_PART)
    {
        return 2*l*l + 2*(m + l);
    }
    else if (part == IMAG_PART)
    {
        return 2*l*l + 2*(m + l) + 1;
    }
}

bool c_set_shc(c_shc const *shc, const PART part, const long l, const long m, const double value)
{
    if (l > shc->bandlimit)
    {
        return false;
    }
    else
    {
        if (part == REAL_PART)
        {
            shc->coefficients[2*l*l + 2*(m + l) ] = value;
            return true;
        }
        else if (part == IMAG_PART)
        {
            shc->coefficients[2*l*l + 2*(m + l) + 1] = value;
            return true;
        }
    }
}


void r_normalize_shc_in_place(r_shc * const shc)
{
    /* 
    Normalizes the shc so that its power spectrum is all ones. 
    It normalizes in place.
    */

    double norm;
    double coeff;

    for (long l = 0; l<=shc->bandlimit; ++l)
    {
        // Calcualte the norm of the coeffiicents of order l
        coeff = r_get_shc(shc, REAL_PART, l, 0);
        norm = coeff*coeff;
        for (long m = 1; m<=l; ++m)
        {
            coeff = r_get_shc(shc, REAL_PART, l, m);
            norm += 2*coeff*coeff;
            
            coeff = r_get_shc(shc, IMAG_PART, l, m);
            norm += 2*coeff*coeff;
        }
        norm = sqrt(norm);

        // Normalize the coefficients
        shc->coefficients[r_lm_to_index(REAL_PART, l, 0)] /= norm;
        for (long m = 1; m<=l; ++m)
        {
            shc->coefficients[r_lm_to_index(REAL_PART, l, m)] /= norm;
            shc->coefficients[r_lm_to_index(IMAG_PART, l, m)] /= norm;
        }
    }
}


void c_normalize_shc_in_place(c_shc * const shc)
{
    /* 
    Normalizes the shc so that its power spectrum is all ones. 
    It normalizes in place.
    */

    double norm;
    double coeff;

    for (long l = 0; l<=shc->bandlimit; ++l)
    {
        // Calculate the norm of the coefficients of order l
        norm = 0;
        for (long m = -l; m<=l; ++m)
        {
            coeff = c_get_shc(shc, REAL_PART, l, m);
            norm += coeff*coeff;

            coeff = c_get_shc(shc, IMAG_PART, l, m);
            norm += coeff*coeff;
        }
        norm = sqrt(norm);

        // Normalize the coefficients
        for (long m = -l; m<=l; ++m)
        {
            shc->coefficients[c_lm_to_index(REAL_PART, l, m)] /= norm;
            shc->coefficients[c_lm_to_index(IMAG_PART, l, m)] /= norm;
        }
    }
}


void r_normalize_shc(r_shc * const shc, r_shc *output_shc)
{
    /* 
    Normalizes the shc so that its power spectrum is all ones. 
    */

    double norm;
    double coeff;

    for (long l = 0; l<=shc->bandlimit; ++l)
    {
        // Calcualte the norm of the coeffiicents of order l
        coeff = r_get_shc(shc, REAL_PART, l, 0);
        norm = coeff*coeff;
        for (long m = 1; m<=l; ++m)
        {
            coeff = r_get_shc(shc, REAL_PART, l, m);
            norm += 2*coeff*coeff;
            
            coeff = r_get_shc(shc, IMAG_PART, l, m);
            norm += 2*coeff*coeff;
        }
        norm = sqrt(norm);

        // Normalize the coefficients
        output_shc->coefficients[r_lm_to_index(REAL_PART, l, 0)] 
            = shc->coefficients[r_lm_to_index(REAL_PART, l, 0)]/norm;
        for (long m = 1; m<=l; ++m)
        {
            output_shc->coefficients[r_lm_to_index(REAL_PART, l, m)] 
                = shc->coefficients[r_lm_to_index(REAL_PART, l, m)]/norm;
            output_shc->coefficients[r_lm_to_index(IMAG_PART, l, m)] 
                = shc->coefficients[r_lm_to_index(IMAG_PART, l, m)]/norm;
        }
    }
}


void c_normalize_shc(c_shc * const shc, c_shc *output_shc)
{
    /* 
    Normalizes the shc so that its power spectrum is all ones. 
    */

    double norm;
    double coeff;

    for (long l = 0; l<=shc->bandlimit; ++l)
    {
        // Calculate the norm of the coefficients of order l
        norm = 0;
        for (long m = -l; m<=l; ++m)
        {
            coeff = c_get_shc(shc, REAL_PART, l, m);
            norm += coeff*coeff;

            coeff = c_get_shc(shc, IMAG_PART, l, m);
            norm += coeff*coeff;
        }
        norm = sqrt(norm);

        // Normalize the coefficients
        for (long m = -l; m<=l; ++m)
        {
            output_shc->coefficients[c_lm_to_index(REAL_PART, l, m)] 
                = shc->coefficients[c_lm_to_index(REAL_PART, l, m)]/norm;
            output_shc->coefficients[c_lm_to_index(IMAG_PART, l, m)] 
                = shc->coefficients[c_lm_to_index(IMAG_PART, l, m)]/norm;
        }
    }
}


double sample_normal(sfmt_t * const sfmt)
{
    /* 
    Generate a double, distributed in the standard normal distribution. 
    REFERENCE:
        https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
    */
    
    double u1 = sfmt_genrand_real3(sfmt);
    double u2 = sfmt_genrand_real3(sfmt);

    return sqrt(-2*log(u1))*cos(2*PI*u2);
}


void r_random_normalized_shc(sfmt_t *sfmt, r_shc *output_shc)
{
    /* 
    REFRENCE:
        http://extremelearning.com.au/how-to-generate-uniformly-random-points-on-n-spheres-and-n-balls/
    */
   
    for (long l = 0; l<=output_shc->bandlimit; ++l)
    {
        r_set_shc(output_shc, REAL_PART, l, 0, sample_normal(sfmt));
        for (long m = 1; m<=l; ++m)
        {
            r_set_shc(output_shc, REAL_PART, l, m, sample_normal(sfmt)/sqrt(2));
            r_set_shc(output_shc, IMAG_PART, l, m, sample_normal(sfmt)/sqrt(2));
        }
    }
    r_normalize_shc_in_place(output_shc);
}


void c_random_normalized_shc(sfmt_t * const sfmt, c_shc *output_shc)
{
    /* 
    REFRENCE:
        http://extremelearning.com.au/how-to-generate-uniformly-random-points-on-n-spheres-and-n-balls/
    */
   
    for (long l = 0; l<=output_shc->bandlimit; ++l)
    {
        // c_set_shc(output_shc, REAL_PART, l, 0, sample_normal(sfmt));
        for (long m = -l; m<=l; ++m)
        {
            c_set_shc(output_shc, REAL_PART, l, m, sample_normal(sfmt));
            c_set_shc(output_shc, IMAG_PART, l, m, sample_normal(sfmt));
        }
    }
    c_normalize_shc_in_place(output_shc);
}


void r_print_shc(r_shc * const shc)
{
    // Print the table header
    printf("Bandlimit = %d.\n", shc->bandlimit);
    printf("=========================================================================\n");
    printf("\t(  l,  m)\t\tReal part\t\tImag part\n");
    printf("=========================================================================\n");

    // Print the coefficients themselves
    for (long l = 0; l<=shc->bandlimit; ++l)
    {
        printf("\t(%3d,%3d)\t\t% 15.15f\n", l, 0, r_get_shc(shc, REAL_PART, l, 0));
        
        for (long m = 1; m<=l; ++m)
        {
            printf("\t(%3d,%3d)\t\t% 15.15f\t% 15.15f\n", 
                    l, m, 
                    r_get_shc(shc, REAL_PART, l, m), 
                    r_get_shc(shc, IMAG_PART, l, m));
        }
        printf("\t----------\n");
    }
}


void c_print_shc(c_shc * const shc)
{
    // Print the table header
    printf("Bandlimit = %d.\n", shc->bandlimit);
    printf("=========================================================================\n");
    printf("\t(  l,  m)\t\tReal part\t\tImag part\n");
    printf("=========================================================================\n");

    // Print the coefficients themselves
    for (long l = 0; l<=shc->bandlimit; ++l)
    {
        for (long m = -l; m<=l; ++m)
        {
            printf("\t(%3d,%3d)\t\t% 15.15f\t% 15.15f\n", 
                    l, m, 
                    c_get_shc(shc, REAL_PART, l, m), 
                    c_get_shc(shc, IMAG_PART, l, m));
        }
        printf("\t----------\n");
    }
}


void r_add_shc(r_shc * restrict shc1, r_shc * restrict shc2, const double alpha)
{
    /* 
    Calculate shc1+alpha*shc2. 
    The code assumes shc1->bandlimit == shc2->bandlimit.
    */

   for (long l = 0; l<shc1->bandlimit; ++l)
   {
       shc1->coefficients[r_lm_to_index(REAL_PART, l, 0)] += alpha*(shc2->coefficients[r_lm_to_index(REAL_PART, l, 0)]);
       for (long m = 1; m<=l; ++m)
       {
           shc1->coefficients[r_lm_to_index(REAL_PART, l, m)] += alpha*(shc2->coefficients[r_lm_to_index(REAL_PART, l, m)]);
           shc1->coefficients[r_lm_to_index(IMAG_PART, l, m)] += alpha*(shc2->coefficients[r_lm_to_index(IMAG_PART, l, m)]);
       }
   }
}


void c_add_shc(c_shc * restrict shc1, c_shc * restrict shc2, const double alpha)
{
    /* 
    Calculate shc1+alpha*shc2. 
    The code assumes shc1->bandlimit == shc2->bandlimit.
    */

   for (long l = 0; l<shc1->bandlimit; ++l)
   {
       for (long m = -l; m<=l; ++m)
       {
           shc1->coefficients[c_lm_to_index(REAL_PART, l, m)] += alpha*(shc2->coefficients[c_lm_to_index(REAL_PART, l, m)]);
           shc1->coefficients[c_lm_to_index(IMAG_PART, l, m)] += alpha*(shc2->coefficients[c_lm_to_index(IMAG_PART, l, m)]);
       }
   }
}


void r_multiply_shc_in_place(r_shc * restrict shc, const double alpha)
{
    /* 
    Calculate alpha*shc. 
    */

   for (long l = 0; l<shc->bandlimit; ++l)
   {
       shc->coefficients[r_lm_to_index(REAL_PART, l, 0)] *= alpha;
       for (long m = 1; m<=l; ++m)
       {
           shc->coefficients[r_lm_to_index(REAL_PART, l, m)] *= alpha;
       }
   }
}


void c_multiply_shc_in_place(c_shc * restrict shc, const double alpha)
{
    /* 
    Calculate shc1+alpha*shc. 
    */

   for (long l = 0; l<shc->bandlimit; ++l)
   {
       for (long m = -l; m<=l; ++m)
       {
           shc->coefficients[c_lm_to_index(REAL_PART, l, m)] *= alpha;
           shc->coefficients[c_lm_to_index(IMAG_PART, l, m)] *= alpha;
       }
   }
}


void c_rotate_spherical_harmonics()
{
    // TODO
}


void r_rotate_spherical_harmonics(r_shc * const restrict shc, double * const restrict rotation, tdesign_cart *td, r_shc * const restrict output_shc)
{
    /* 
    Calcualte the spherical harmonics coefficients of the bandlimited function with spherical harmonics coefficients shc.
    rotation is a an array of size 4, which is a quaternion representation of a rotation.
    Result is saved in output_shc, which is assumed to be all zeros..

    Must run alegendre_init once before using this function.
    This is a relatively naive implementation, without many optimizations.
    */

    double theta, phi;
    double rotated_x[3];

    double sh_rp, sh_ip;
    double func_val;

    const double integration_factor = 2*PI/td->length;

    for (size_t i = 0; i<td->length; ++i)
    {
        apply_rotation(&(td->tdesign[3*i]), rotation, rotated_x);
        cartesian_unit_vector_to_spherical(rotated_x, &theta, &phi);
        
        r_eval_sf(shc, theta, phi, &func_val);

        for (long l = 0; l<shc->bandlimit; ++l)
        {
            for (long m = 1; m<=l; ++m)
            {
                eval_sh(l, m, theta, phi, &sh_rp, &sh_ip);
                output_shc->coefficients[r_lm_to_index(REAL_PART, l, m)] += integration_factor*(m%2==0 ? 1 : -1)*func_val*sh_rp;
                output_shc->coefficients[r_lm_to_index(IMAG_PART, l, m)] += integration_factor*(m%2==0 ? 1 : -1)*func_val*sh_ip;
            }
            // Handle m =0

            eval_sh(l, 0, theta, phi, &sh_rp, &sh_ip);
            output_shc->coefficients[r_lm_to_index(REAL_PART, l, 0)] += integration_factor*func_val*sh_rp;
            output_shc->coefficients[r_lm_to_index(IMAG_PART, l, 0)] += integration_factor*func_val*sh_ip;
        }
    }
}

void r_eval_sf(r_shc * const shc, const double theta, const double phi, double * const restrict real_part)
{
    /* 
    Evaluate real-valued spherical function represented by shc.
    
    Must run alegendre_init once before using this function.
    */

   double sh_rp, sh_ip;
   *real_part = 0.0;
   for (long l = 0; l<=shc->bandlimit; ++l)
   {
       for (long m = 1; m<=l; ++m)
       {
           eval_sh(l, m, theta, phi, &sh_rp, &sh_ip);
           *real_part   += r_get_shc(shc, REAL_PART, l, m)*sh_rp 
                            - r_get_shc(shc, IMAG_PART, l, m)*sh_ip;
       }
       *real_part *= 2;

        eval_sh(l, 0, theta, phi, &sh_rp, &sh_ip);
       *real_part += r_get_shc(shc, REAL_PART, l, 0)*sh_rp;
   }
}

double c_eval_sf(c_shc * const shc, const double theta, const double phi, const PART part, double * const restrict real_part, double * const restrict imag_part)
{
    // Evaluate complex-valued spherical function represented by shc
}

void eval_sh(const long l, const long m, const double theta, const double phi, double * const restrict real_part, double * const restrict imag_part)
{
    /* 
    Evaluate the real part or the imaginary part of the spherical harmonics 
        Y_{l}^{m}(theta,phi) := e^{i m phi} * S_{l}^{m} (theta)
    where S_{l}^{m} was defined in alegendre() docs in alegendre.h.

    Output is placed in real_part and imag_part.

    Must run alegendre_init once before using this function.
    */

    *real_part = alegendre(l, m, theta);
    *imag_part = *real_part;
    
    *real_part *= cos(m*phi);
    *imag_part *= sin(m*phi);
}

void cartesian_unit_vector_to_spherical(double * const restrict x, double * const restrict theta, double * const restrict phi)
{
    /* 
    x is an array of length 3, representing a unit vector in R^3.
    Converts x to spherical coordinates (theta, phi) where 
        0<=theta<=pi and 0<=phi<2pi.

    Code performs no input checks.
    */
   *phi = atan2(x[1], x[0]);
   if (*phi<0)
   {
       *phi += 2*PI;
   }
   *theta = acos(x[2]);
}