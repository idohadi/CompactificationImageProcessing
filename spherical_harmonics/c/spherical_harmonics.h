// TODO

#include <stdint.h>
#include <stdbool.h>
#include "SFMT.h"

/* Typedefs and enum's */

/*  
    for real-valued spherical functions: spherical harmonics coefficients format is
        f_{r, 0, 0}, 
            f_{r, 1, 0}, f_{r, 1, 1}, f_{i, 1, 1}, 
                f_{r, 2, 0}, f_{r, 2, 1}, f_{i, 2, 1}, f_{r, 2, 2}, f_{i, 2, 2}, ...
                    f_{r, bandlimit, 0}, f_{r, bandlimit, 1}, f_{i, bandlimit, 1}, ..., f_{r, bandlimit, bandlimit}, f_{i, bandlimit, bandlimit}
*/
typedef 
    struct r_shc
    {
        double *coefficients;
        size_t bandlimit;
    }
    r_shc;

/*  
    for complex-valued spherical functions: spherical harmonics coefficients format is
        f_{r, 0, 0}, 
            f_{r, 1, -1}, f_{i, 1, -1}, f_{r, 1, 0}, f_{i, 1, 0}, f_{r, 1, 1}, f_{i, 1, 1}, 
                f_{r, 2, -2}, f_{i, 2, -2}, f_{r, 2, -1}, f_{i, 2, -1}, f_{r, 2, 0}, f_{i, 2, 0}, f_{r, 2, 1}, f_{i, 2, 1}, f_{r, 2, 2}, f_{i, 2, 2}, 
                ...
*/
typedef 
    struct c_shc
    {
        double *coefficients;
        size_t bandlimit;
    }
    c_shc;

typedef enum PART {REAL_PART, IMAG_PART} PART;


/* Functions */
double *r_allocate_coefficients(const size_t bandlimit);


double *c_allocate_coefficients(const size_t bandlimit);


void r_deallocate_coefficients(double * const coefficients);


void c_deallocate_coefficients(double * const coefficients);


r_shc r_init_shc(const size_t bandlimit, double * const coefficients);


c_shc c_init_shc(const size_t bandlimit, double * const coefficients);


r_shc r_init_shc2(const size_t bandlimit);


c_shc c_init_shc2(const size_t bandlimit);


double r_get_shc(r_shc const *shc, const PART part, const long l, const long m);


double c_get_shc(c_shc const *shc, const PART part, const long l, const long m);


bool r_set_shc(r_shc const *shc, const PART part, const long l, const long m, const double value);


long r_lm_to_index(const PART part, const long l, const long m);


long c_lm_to_index(const PART part, const long l, const long m);


bool c_set_shc(c_shc const *shc, const PART part, const long l, const long m, const double value);


void r_normalize_shc_in_place(r_shc * const shc);


void c_normalize_shc_in_place(c_shc * const shc);


void r_normalize_shc(r_shc * const shc, r_shc *output_shc);


void c_normalize_shc(c_shc * const shc, c_shc *output_shc);


double sample_normal(sfmt_t * const sfmt);


void r_random_normalized_shc(sfmt_t *sfmt, r_shc *output_shc);


void c_random_normalized_shc(sfmt_t * const sfmt, r_shc * const output_shc);


void r_print_shc(r_shc * const shc);


void c_print_shc(c_shc * const shc);


void r_add_shc(r_shc * restrict shc1, r_shc * restrict shc2, const double alpha);


void c_add_shc(c_shc * restrict shc1, c_shc * restrict shc2, const double alpha);


void r_multiply_shc_in_place(r_shc * restrict shc, const double alpha);


void c_multiply_shc_in_place(c_shc * restrict shc, const double alpha);