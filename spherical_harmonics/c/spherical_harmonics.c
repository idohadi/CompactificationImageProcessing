// TODO

#include <stdbool.h>
#include <math.h>
#include "spherical_harmonics.h"

double *r_allocate_coefficients()
{
    // TODO
    // Allocate memory for SHC
}

double *c_allocate_coefficients()
{
    // TODO
    // Allocate memory for SHC
}

double *r_deallocate_coefficients()
{
    // TODO
}

double *c_deallocate_coefficients()
{
    // TODO
}

void r_init_shc()
{
    // TODO: initialize struct containing the bandlimit and coefficients
    // coefficients are given as pointer to array of appropriate sizeZ
}

void c_init_shc()
{
    // TODO: initialize struct containing the bandlimit and coefficients
    // coefficients are given as pointer to array of appropriate sizeZ
}

void r_init_shc2()
{
    // TODO: initialize struct containing the bandlimit and coefficients
    // space for coefficients is allocated and they are initialized to all zeros
}

void c_init_shc2()
{
    // TODO: initialize struct containing the bandlimit and coefficients
    // space for coefficients is allocated and they are initialized to all zeros
}

void r_destroy_shc()
{
    // TODO: destroys struct containing the bandlimit and coefficients
    // doesn't deallocate coefficients memory
}

void c_destroy_shc()
{
    // TODO: initialize struct containing the bandlimit and coefficients
    // doesn't deallocate coefficients memory
}

void r_destroy_shc2()
{
    // TODO: initialize struct containing the bandlimit and coefficients
    // deallocates coefficients memory
}

void c_destroy_shc2()
{
    // TODO: initialize struct containing the bandlimit and coefficients
    // deallocates coefficients memory
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


void r_normalize_shc(const r_shc *shc)
{
    // It normalizes in place
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


void c_normalize_shc(const c_shc *shc)
{
    // It normalizes in place
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


void r_random_normalized_shc()
{
    // TODO
}


void c_random_normalized_shc()
{
    // TODO
}


void r_print_shc()
{
    // TODO
}


void c_print_shc()
{
    // TODO
}


void r_add_shc()
{
    // TODO: calcualte shc1+alpha*shc2
    // Should use restricted pointers
}


void c_add_shc()
{
    // TODO: calcualte shc1+alpha*shc2
    // Should use restricted pointers
}


void r_multiply_shc_in_place()
{
    // TODO: calculate alpha*shc2
}


void c_multiply_shc_in_place()
{
    // TODO: calculate alpha*shc2
}


void c_rotate_spherical_harmonics()
{
    // TODO
}


void r_rotate_spherical_harmonics()
{
    // TODO
}

