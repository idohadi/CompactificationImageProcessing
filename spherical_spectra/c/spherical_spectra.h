#include <stdint.h>
#include "clebsch_gordan_coefficients.h"
#include "spherical_harmonics.h"


typedef size_t *** r_bispectrum_lookout_table;
typedef r_bispectrum_lookout_table r_blt;

typedef size_t **** c_bispectrum_lookout_table;
typedef c_bispectrum_lookout_table c_blt;


c_blt c_build_bispectrum_lookup_table(size_t bandlimit);


void c_destroy_bispectrum_lookup_table(c_blt table, size_t bandlimit);


r_bispectrum_lookout_table r_build_bispectrum_lookup_table(size_t bandlimit);


void r_destroy_bispectrum_lookup_table(r_bispectrum_lookout_table table, size_t bandlimit);


/* Power spectrum fucntions */
double *c_allocate_power_spectrum(const size_t bandlimit);


double *r_allocate_power_spectrum(const size_t bandlimit);


void c_power_spectrum(c_shc * const shc, double *c_power_spectrum);


void r_power_spectrum(r_shc * const shc, double *r_pow_spec);


/* Bispectrum functions */

void c_bispectrum(c_shc * const shc, const c_blt lookup, cg_table * const table, double *c_bisp);


double r_bispectral_invariant_real_part(r_shc * const shc, const long l1, const long l2, const long l, cg_table * const table);


double r_bispectral_invariant_imaginary_part(r_shc * const shc, const long l1, const long l2, const long l, cg_table * const table);


double c_bispectral_invariant_real_part(c_shc * const shc, const long l1, const long l2, const long l, const cg_table *table);


double c_bispectral_invariant_imaginary_part(c_shc * const shc, const long l1, const long l2, const long l, const cg_table *table);


void r_bispectrum(r_shc * const shc, const c_blt lookup, cg_table * const table, double *r_bisp);


void r_bispectrum_gradient(r_shc * const shc, const c_blt lookup, cg_table * const table, double *r_bisp_grad);


void c_bispectrum_gradient(c_shc * const shc, const c_blt lookup, cg_table * const table, double *c_bisp_grad);


/* Printing functions */
void r_print_power_spectrum(double *r_pow_spec, size_t bandlimit);


void c_print_power_spectrum(double *c_pow_spec, size_t bandlimit);


void r_print_all_bispectral_invariants(r_shc * const shc, cg_table * const table);