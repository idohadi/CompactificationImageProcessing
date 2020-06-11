// TODO: write docs in doxygen

#include <stdbool.h>
#include <math.h>
#include <stdio.h>
#include "spherical_harmonics.h"
#include "SFMT.h"

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
            r_set_shc(&shc, REAL_PART, l, m, 0.0);
            r_set_shc(&shc, IMAG_PART, l, m, 0.0);
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
            shc->coefficients[c_lm_to_index(REAL_PART, l, m)] /= norm;
            shc->coefficients[c_lm_to_index(IMAG_PART, l, m)] /= norm;
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
    r_normalize_shc(output_shc);
}


void c_random_normalized_shc(sfmt_t * const sfmt, r_shc * const output_shc)
{
    /* 
    REFRENCE:
        http://extremelearning.com.au/how-to-generate-uniformly-random-points-on-n-spheres-and-n-balls/
    */
   
    for (long l = 0; l<=output_shc->bandlimit; ++l)
    {
        r_set_shc(output_shc, REAL_PART, l, 0, sample_normal(sfmt));
        for (long m = -l; m<=l; ++m)
        {
            r_set_shc(output_shc, REAL_PART, l, m, sample_normal(sfmt));
            r_set_shc(output_shc, IMAG_PART, l, m, sample_normal(sfmt));
        }
    }
    c_normalize_shc(output_shc);
}


void r_print_shc(r_shc * const shc)
{
    // Print the table header
    printf("Bandlimit = %d.", shc->bandlimit);
    printf("===========================================================================================================\n");
    printf("\t(part,order,degree)\t\tValue\n");
    printf("===========================================================================================================\n");

    // Print the coefficients themselves
    for (long l = 0; l<shc->bandlimit; ++l)
    {
        printf("\t(r,%2-d,%-2d)\t\t%-15.15f\n", l, 0, r_get_shc(shc, REAL_PART, l, 0));
        for (long m = 1; m<=l; ++m)
        {
            printf("\t(r,%2-d,%-2d)\t\t%-15.15f\n", l, m, r_get_shc(shc, REAL_PART, l, m));
            printf("\t(i,%2-d,%-2d)\t\t%-15.15f\n", l, m, r_get_shc(shc, IMAG_PART, l, m));
        }
    }
}


void c_print_shc(c_shc * const shc)
{
    // Print the table header
    printf("Bandlimit = %d.", shc->bandlimit);
    printf("===========================================================================================================\n");
    printf("\t(part,order,degree)\t\tValue\n");
    printf("===========================================================================================================\n");

    // Print the coefficients themselves
    for (long l = 0; l<shc->bandlimit; ++l)
    {
        for (long m = -l; m<=l; ++m)
        {
            printf("\t(r,%2-d,%-2d)\t\t%-15.15f\n", l, m, r_get_shc(shc, REAL_PART, l, m));
            printf("\t(i,%2-d,%-2d)\t\t%-15.15f\n", l, m, r_get_shc(shc, IMAG_PART, l, m));
        }
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


void r_rotate_spherical_harmonics()
{
    // TODO
}

