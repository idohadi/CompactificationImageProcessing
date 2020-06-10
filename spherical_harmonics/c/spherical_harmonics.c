// TODO

#include "spherical_harmonics.h"

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

double r_get_shc(r_shc *shc, PART part, long int l, long int m)
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

double c_get_shc()
{
    // TODO
}

void r_set_shc()
{
    // TODO
}

void c_set_shc()
{
    // TODO
}

void r_normalize_shc()
{
    // TODO
}

void c_normalize_shc()
{
    // TODO
}

void r_random_normalized_shc()
{
    // TODO
}

void c_random_normalized_shc()
{
    // TODO
}

double *r_allocate_shc()
{
    // TODO
    // Allocate memory for SHC
}

double *c_allocate_shc()
{
    // TODO
    // Allocate memory for SHC
}

double *r_deallocate_shc()
{
    // TODO
}

double *c_deallocate_shc()
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

