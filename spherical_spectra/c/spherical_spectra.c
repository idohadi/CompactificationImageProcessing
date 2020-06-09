// TODO: write docs
/**
 * prefix c_ means this is a function for spherical harmonics expansion of a complex-valued spherical function
 * prefix r_ means this is a function for spherical harmonics expansion of a real-valued spherical function
 */


#include <stdint.h>
#include <stdio.h>
#include "clebsch_gordan_coefficients.h"
#include "spherical_spectra.h"

c_bispectrum_lookout_table c_build_bispectrum_lookup_table(size_t bandlimit)
{
    // TODO
}

void c_destroy_bispectrum_lookup_table(c_bispectrum_lookout_table table, size_t L)
{
    // TODO
}

r_bispectrum_lookout_table r_build_bispectrum_lookup_table(size_t bandlimit)
{
    size_t index = 0; 

    size_t ***lookup_table = malloc((bandlimit+1)*sizeof(size_t **));
    for (long int l1 = 0; l1<=bandlimit; ++l1)
    {
        lookup_table[l1] = malloc((l1+1)*sizeof(size_t *));
        for (long int l2=0; l2<=l1; ++l2)
        {
            lookup_table[l1][l2] = malloc((minimum(l1+l2, bandlimit) - absolute_value(l1-l2) + 1)*sizeof(size_t));
            for (long int l = absolute_value(l1-l2); l<=minimum(l1+l2, bandlimit); ++l)
            {
                lookup_table[l1][l2][l - absolute_value(l1-l2)] = index;
                ++index;
            }
        }
    }
    return lookup_table;
}

void r_destroy_bispectrum_lookup_table(r_bispectrum_lookout_table table, size_t bandlimit)
{
    for (long int l1 = 0; l1<=bandlimit; ++l1)
    {
        for (long int l2=0; l2<=l1; ++l2)
        {
            free(table[l1][l2]);
        }
        free(table[l1]);
    }
    free(table);
}

double r_get_spherical_harmonics(const double *spherical_harmonics_coeffs, char real_part_or_imaginary_part, long int l, long int m)
{
    /*  
    real_part_or_imaginary_part = 'r' or 'i'
    spherical_harmonics_coeffs format is
        f_{r, 0, 0}, 
            f_{r, 1, 0}, f_{r, 1, 1}, f_{i, 1, 1}, 
                f_{r, 2, 0}, f_{r, 2, 1}, f_{i, 2, 1}, f_{r, 2, 2}, f_{i, 2, 2}, ...
                    f_{r, bandlimit, 0}, f_{r, bandlimit, 1}, f_{i, bandlimit, 1}, ..., f_{r, bandlimit, bandlimit}, f_{i, bandlimit, bandlimit}
    */

   if (m==0)
   {
       if (real_part_or_imaginary_part=='i')
       {
           return 0.0;
       }
       else if (real_part_or_imaginary_part=='r')
       {
           return spherical_harmonics_coeffs[l*l];
       }
   }

   if (m>0)
   {
       if (real_part_or_imaginary_part=='r')
       {
           return spherical_harmonics_coeffs[l*l + 2*m - 1];
       }
       else if (real_part_or_imaginary_part=='i')
       {
           return spherical_harmonics_coeffs[l*l + 2*m];
       }
   }

   if (m<0)
   {
       if (real_part_or_imaginary_part=='r')
       {
           return minus_one_power(m)*spherical_harmonics_coeffs[l*l - 2*m - 1];
       }
       else if (real_part_or_imaginary_part=='i')
       {
           return minus_one_power(m+1)*spherical_harmonics_coeffs[l*l - 2*m];
       }
   }
}

void r_set_spherical_harmonics(double *spherical_harmonics_coeffs, long int l, long int m, double value)
{
    // TODO
}

void c_calculate_bispectrum(const size_t bandlimit, const double *c_spherical_harmonics_coeffs, double *c_bipsectrum)
{
    // TODO
}

double r_calculate_bispectral_invariant_real_part(  double * const r_spherical_harmonics_coeffs, 
                                                    const size_t bandlimit, 
                                                    const long int l1, 
                                                    const long int l2, 
                                                    const long int l, 
                                                    const cg_table clebsch_gordan_coeffs)
{
    double invariant = 0;

    double U = 0;
    double M = 0;

    long int lower_bound = 0;
    long int upper_bound = 0;
    for (long int m = -l; m<=l; ++m)
    {
        U = 0;
        M = 0;

        lower_bound = (-l1)>=(m-l2) ? (-l1) : (m-l2);
        upper_bound = l1>=(m+l2) ? (m+l2) : (l1);
        for (long int m1 = lower_bound; m1<=upper_bound; ++m1)
        {
            U +=    get_clebsch_gordan_coefficient(clebsch_gordan_coeffs, l1, l2, l, m, m1)
                    *( 
                        r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l1, m1)
                        *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l2, m-m1)
                        -   r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l1, m1)
                            *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l2, m-m1)
                    );
            M +=    get_clebsch_gordan_coefficient(clebsch_gordan_coeffs, l1, l2, l, m, m1)
                    *(
                        r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l1, m1)
                        *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l2, m-m1)
                        +   r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l1, m1)
                            *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l2, m-m1)
                    );
        }

        invariant 
            +=  r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l, m)*U 
                +r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l, m)*M;
    }

    return invariant;
}

double r_calculate_bispectral_invariant_imaginary_part( double * const r_spherical_harmonics_coeffs, 
                                                        const size_t bandlimit, 
                                                        const long int l1, 
                                                        const long int l2, 
                                                        const long int l, 
                                                        const cg_table clebsch_gordan_coeffs)
{
    double invariant = 0;

    double U = 0;
    double M = 0;

    long int lower_bound = 0;
    long int upper_bound = 0;
    for (long int m = -l; m<=l; ++m)
    {
        U = 0;
        M = 0;

        lower_bound = (-l1)>=(m-l2) ? (-l1) : (m-l2);
        upper_bound = l1>=(m+l2) ? (m+l2) : (l1);
        for (long int m1 = lower_bound; m1<=upper_bound; ++m1)
        {
            U +=    get_clebsch_gordan_coefficient(clebsch_gordan_coeffs, l1, l2, l, m, m1)
                    *( 
                        r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l1, m1)
                        *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l2, m-m1)
                        -   r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l1, m1)
                            *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l2, m-m1)
                    );
            M +=    get_clebsch_gordan_coefficient(clebsch_gordan_coeffs, l1, l2, l, m, m1)
                    *(
                        r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l1, m1)
                        *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l2, m-m1)
                        +   r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l1, m1)
                            *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l2, m-m1)
                    );
        }

        invariant 
            +=  r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l, m)*U 
                -r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l, m)*M;
    }

    return invariant;
}

void r_calculate_bispectrum(double * const r_spherical_harmonics_coeffs, 
                            const size_t bandlimit, 
                            size_t *** const bispectrum_lookup_table, 
                            const cg_table clebsch_gordan_coeffs, 
                            double *r_bispectrum)
{
    // #pragma omp parallel for collapse(3)
    
    for (long int l1 = 0; l1<=bandlimit; ++l1)
    {
        for (long int l2 = 0; l2<=l1; ++l2)
        {
            for (long int l = l1-l2; l<=l1+l2 && l<=bandlimit; ++l)
            {
                size_t row_index = bispectrum_lookup_table[l1][l2][l-(l1-l2)];

                r_bispectrum[row_index]
                    = r_calculate_bispectral_invariant_real_part(r_spherical_harmonics_coeffs, bandlimit, l1, l2, l, clebsch_gordan_coeffs);
            }
        }
    }
}

void c_calculate_bispectrum_gradient(   const double *c_spherical_harmonics_coeffs, 
                                        const size_t bandlimit, 
                                        const size_t ***bispectrum_lookup_table, 
                                        const cg_table clebsch_gordan_coeffs, 
                                        double *c_bipsectrum_gradient)
{
    // TODO
}

void r_calculate_bispectrum_gradient(   double * const r_spherical_harmonics_coeffs, 
                                        const size_t bandlimit, 
                                        size_t *** const bispectrum_lookup_table, 
                                        const cg_table clebsch_gordan_coeffs, 
                                        double *r_bipsectrum_gradient)
{
    /* I assume the array r_bipsectrum_gradient is initialized with zeros.
    TODO: document the hell out of this, including an external file explaining the derivation. */
    for (long int l1 = 0; l1<=bandlimit; ++l1)
    {
        for (long int l2 = 0; l2<=l1; ++l2)
        {
            for (long int l = l1-l2; l<=l1+l2 && l<=bandlimit; ++l)
            {
                double Uplus = 0;
                double Uminus = 0;
                double Mplus = 0;
                double Mminus = 0;

                // Working on the derivative of b_{l1,l2,l}
                size_t row_index = bispectrum_lookup_table[l1][l2][l-(l1-l2)]*(bandlimit+1)*(bandlimit+1);

                // *********************************
                // Derivative w.r.t f_{sigma,l,m}
                // *********************************
                long int lower_bound = 0;
                long int upper_bound = 0;
                
                for (long int m = 1; m<=l; ++m) // for m>=1
                {
                    Uplus = 0;
                    Uminus = 0;
                    Mplus = 0;
                    Mminus = 0;

                    lower_bound = (-l1)>=(m-l2) ? (-l1) : (m-l2);
                    upper_bound = l1>=(m+l2) ? (m+l2) : (l1);
                    for (long int m1 = lower_bound; m1<=upper_bound; ++m1)
                    {
                        Uplus +=    get_clebsch_gordan_coefficient(clebsch_gordan_coeffs, l1, l2, l, m, m1)
                                *( 
                                    r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l1, m1)
                                    *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l2, m-m1)
                                    -   r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l1, m1)
                                        *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l2, m-m1)
                                );
                        Mplus +=    get_clebsch_gordan_coefficient(clebsch_gordan_coeffs, l1, l2, l, m, m1)
                                *(
                                    r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l1, m1)
                                    *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l2, m-m1)
                                    +   r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l1, m1)
                                        *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l2, m-m1)
                                );
                    }

                    lower_bound = (-l1)>=(-m-l2) ? (-l1) : (-m-l2);
                    upper_bound = l1>=(-m+l2) ? (-m+l2) : (l1);
                    for (long int m1 = lower_bound; m1<=upper_bound; ++m1)
                    {
                        Uminus +=    get_clebsch_gordan_coefficient(clebsch_gordan_coeffs, l1, l2, l, -m, m1)
                                *( 
                                    r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l1, m1)
                                    *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l2, -m-m1)
                                    -   r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l1, m1)
                                        *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l2, -m-m1)
                                );
                        Mminus +=    get_clebsch_gordan_coefficient(clebsch_gordan_coeffs, l1, l2, l, -m, m1)
                                *(
                                    r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l1, m1)
                                    *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l2, -m-m1)
                                    +   r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l1, m1)
                                        *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l2, -m-m1)
                                );
                    }

                    // The derivative w.r.t f_{r,l,m} (m>=1)
                    r_bipsectrum_gradient[row_index+l*l+2*m-1] += Uplus + ((m % 2)==0 ? 1 : -1)*Uminus;

                    // The derivative w.r.t f_{i,l,m} (m>=1)
                    r_bipsectrum_gradient[row_index+l*l+2*m] += Mplus - ((m % 2)==0 ? 1 : -1)*Mminus;
                }

                // The derivative w.r.t f_{r,l,0}
                lower_bound = (-l1)>=(-l2) ? (-l1) : (-l2);
                upper_bound = l1>=l2 ? (l2) : (l1);
                for (long int m1 = lower_bound; m1<=upper_bound; ++m1)
                {
                    r_bipsectrum_gradient[row_index+l*l] +=    get_clebsch_gordan_coefficient(clebsch_gordan_coeffs, l1, l2, l, 0, m1)
                                            *( 
                                                r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l1, m1)
                                                *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l2, -m1)
                                                -   r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l1, m1)
                                                    *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l2, -m1)
                                            );
                }
                

                // *********************************
                // Derivative w.r.t f_{sigma,l1,m1}
                // *********************************
                lower_bound = 0;
                upper_bound = 0;
                
                for (long int m1 = 1; m1<=l1; ++m1) // for m1>=1
                {
                    Uplus = 0;
                    Uminus = 0;
                    Mplus = 0;
                    Mminus = 0;

                    lower_bound = (-l)>=(m1-l2) ? (-l) : (m1-l2);
                    upper_bound = l>=(m1+l2) ? (m1+l2) : l;
                    for (long int m = lower_bound; m<=upper_bound; ++m)
                    {
                        Uplus +=    get_clebsch_gordan_coefficient(clebsch_gordan_coeffs, l1, l2, l, m, m1)
                                *( 
                                    r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l, m)
                                    *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l2, m-m1)
                                    +   r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l, m)
                                        *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l2, m-m1)
                                );
                        Mplus +=    get_clebsch_gordan_coefficient(clebsch_gordan_coeffs, l1, l2, l, m, m1)
                                *(
                                    r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l, m)
                                    *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l2, m-m1)
                                    -   r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l, m)
                                        *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l2, m-m1)
                                );
                    }

                    lower_bound = (-l)>=(-m1-l2) ? (-l) : (-m1-l2);
                    upper_bound = l>=(-m1+l2) ? (-m1+l2) : l;
                    for (long int m = lower_bound; m<=upper_bound; ++m)
                    {
                        Uminus +=    get_clebsch_gordan_coefficient(clebsch_gordan_coeffs, l1, l2, l, m, -m1)
                                *( 
                                    r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l, m)
                                    *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l2, m+m1)
                                    +   r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l, m)
                                        *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l2, m+m1)
                                );
                        Mminus +=    get_clebsch_gordan_coefficient(clebsch_gordan_coeffs, l1, l2, l, m, -m1)
                                *(
                                    r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l, m)
                                    *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l2, m+m1)
                                    -   r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l, m)
                                        *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l2, m+m1)
                                );
                    }

                    // The derivative w.r.t f_{r,l1,m1} (m1>=1)
                    r_bipsectrum_gradient[row_index+l1*l1+2*m1-1] += Uplus + ((m1 % 2)==0 ? 1 : -1)*Uminus;
                
                    // The derivative w.r.t f_{i,l1,m1} (m1>=1)
                    r_bipsectrum_gradient[row_index+l1*l1+2*m1] += Mplus - ((m1 % 2)==0 ? 1 : -1)*Mminus;
                }

                // The derivative w.r.t f_{r,l1,0}
                lower_bound = (-l)>=(-l2) ? (-l) : (-l2);
                upper_bound = l>=l2 ? l2 : l;
                for (long int m = lower_bound; m<=upper_bound; ++m)
                {
                    r_bipsectrum_gradient[row_index+l1*l1] +=    get_clebsch_gordan_coefficient(clebsch_gordan_coeffs, l1, l2, l, m, 0)
                            *( 
                                r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l, m)
                                *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l2, m)
                                +   r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l, m)
                                    *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l2, m)
                            );
                }
                


                // *********************************
                // Derivative w.r.t f_{sigma,l2,m2}
                // *********************************
                lower_bound = 0;
                upper_bound = 0;
                
                for (long int m2 = 1; m2<=l2; ++m2) // for m2>=1
                {
                    Uplus = 0;
                    Uminus = 0;
                    Mplus = 0;
                    Mminus = 0;

                    lower_bound = (-l)>=(m2-l1) ? (-l) : (m2-l1);
                    upper_bound = l>=(m2+l1) ? (m2+l1) : l;
                    for (long int m = lower_bound; m<=upper_bound; ++m)
                    {
                        Uplus +=    get_clebsch_gordan_coefficient(clebsch_gordan_coeffs, l1, l2, l, m, m-m2)
                                *( 
                                    r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l, m)
                                    *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l1, m-m2)
                                    +   r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l, m)
                                        *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l1, m-m2)
                                );
                        Mplus +=    get_clebsch_gordan_coefficient(clebsch_gordan_coeffs, l1, l2, l, m, m-m2)
                                *(
                                    r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l, m)
                                    *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l1, m-m2)
                                    -   r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l, m)
                                        *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l1, m-m2)
                                );
                    }

                    lower_bound = (-l)>=(-m2-l1) ? (-l) : (-m2-l1);
                    upper_bound = l>=(-m2+l1) ? (-m2+l1) : l;
                    for (long int m = lower_bound; m<=upper_bound; ++m)
                    {
                        Uminus +=    get_clebsch_gordan_coefficient(clebsch_gordan_coeffs, l1, l2, l, m, m+m2)
                                *( 
                                    r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l, m)
                                    *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l1, m+m2)
                                    +   r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l, m)
                                        *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l1, m+m2)
                                );
                        Mminus +=    get_clebsch_gordan_coefficient(clebsch_gordan_coeffs, l1, l2, l, m, m+m2)
                                *(
                                    r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l, m)
                                    *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l1, m+m2)
                                    -   r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l, m)
                                        *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l1, m+m2)
                                );
                    }

                    // The derivative w.r.t f_{r,l2,m2} (m2>=1)
                    r_bipsectrum_gradient[row_index+l2*l2+2*m2-1] += Uplus + ((m2 % 2)==0 ? 1 : -1)*Uminus;
                
                    // The derivative w.r.t f_{i,l2,m2} (m2>=1)
                    r_bipsectrum_gradient[row_index+l2*l2+2*m2] += Mplus - ((m2 % 2)==0 ? 1 : -1)*Mminus;
                }

                // The derivative w.r.t f_{r,l2,0}
                lower_bound = (-l)>=(-l1) ? (-l) : (-l1);
                upper_bound = l>=l1 ? l1 : l;
                for (long int m = lower_bound; m<=upper_bound; ++m)
                {
                    r_bipsectrum_gradient[row_index+l2*l2] +=    get_clebsch_gordan_coefficient(clebsch_gordan_coeffs, l1, l2, l, m, m)
                            *( 
                                r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l, m)
                                *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'r', l1, m)
                                +   r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l, m)
                                    *r_get_spherical_harmonics(r_spherical_harmonics_coeffs, 'i', l1, m)
                            );
                }
            }
        }
    }
}

void c_calculate_power_spectrum(const size_t bandlimit, const double *c_spherical_harmonics_coeffs, double *r_power_spectrum)
{
    // TODO
}

void r_calculate_power_spectrum(const size_t bandlimit, const double *r_spherical_harmonics_coeffs, double *c_power_spectrum)
{
    // TODO
}
