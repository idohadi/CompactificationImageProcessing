// TODO
#include "clebsch_gordan_coefficients.h"
#include <stdint.h>

long int maximum(long int a, long int b);

long int minimum(long int a, long int b);

long int absolute_value(long int a);

double minus_one_power(long int m);

size_t ***r_build_bispectrum_lookup_table(size_t bandlimit);

void r_destroy_bispectrum_lookup_table(size_t ***lookup_table, size_t bandlimit);

double r_get_spherical_harmonics(const double *spherical_harmonics_coeffs, char real_part_or_imaginary_part, long int l, long int m);

double unmixed_sum_utility_function(const double *r_spherical_harmonics_coeffs, 
                                    long int l1, 
                                    long int l2, 
                                    long int l, 
                                    long int m, 
                                    const cb_table clebsch_gordan_coeffs);

double mixed_sum_utility_function(  const double *r_spherical_harmonics_coeffs, 
                                    long int l1, 
                                    long int l2, 
                                    long int l, 
                                    long int m, 
                                    const cb_table clebsch_gordan_coeffs);

void r_calculate_bispectrum(double * const r_spherical_harmonics_coeffs, 
                            const size_t bandlimit, 
                            size_t *** const bispectrum_lookup_table, 
                            const cb_table clebsch_gordan_coeffs, 
                            double *r_bispectrum);

void r_calculate_bispectrum_gradient(   double * const r_spherical_harmonics_coeffs, 
                                        const size_t bandlimit, 
                                        size_t *** const bispectrum_lookup_table, 
                                        const cb_table clebsch_gordan_coeffs, 
                                        double *r_bipsectrum_gradient);

double r_calculate_bispectral_invariant_imaginary_part( double * const r_spherical_harmonics_coeffs, 
                                                        const size_t bandlimit, 
                                                        const long int l1, 
                                                        const long int l2, 
                                                        const long int l, 
                                                        const cb_table clebsch_gordan_coeffs);

double r_calculate_bispectral_invariant_real_part(  double * const r_spherical_harmonics_coeffs, 
                                                    const size_t bandlimit, 
                                                    const long int l1, 
                                                    const long int l2, 
                                                    const long int l, 
                                                    const cb_table clebsch_gordan_coeffs);
