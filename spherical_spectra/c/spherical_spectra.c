// TODO: write docs in doxygen
/**
 * prefix c_ means this is a function for spherical harmonics expansion of a complex-valued spherical function
 * prefix r_ means this is a function for spherical harmonics expansion of a real-valued spherical function
 */


#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "spherical_spectra.h"


/* Lookup table functions */
c_blt c_build_bispectrum_lookup_table(size_t bandlimit)
{
    size_t index = 0; 

    size_t ****lookup_table = malloc((bandlimit+1)*sizeof(size_t ***));
    for (long l1 = 0; l1<=bandlimit; ++l1)
    {
        lookup_table[l1] = malloc((l1+1)*sizeof(size_t **));
        for (long l2=0; l2<=l1; ++l2)
        {
            lookup_table[l1][l2] = malloc((((l1+l2)>=bandlimit ? bandlimit : l1+l2)  - (l1>=l2 ? l1-l2 : l2-l1) + 1)*sizeof(size_t *));
            for (long l = (l1>=l2 ? l1-l2 : l2-l1); l<=l1+l2 && l<=bandlimit; ++l)
            {
                lookup_table[l1][l2][l - (l1>=l2 ? l1-l2 : l2-l1)] = malloc(2*sizeof(size_t));
                lookup_table[l1][l2][l - (l1>=l2 ? l1-l2 : l2-l1)][0] = index++;
                lookup_table[l1][l2][l - (l1>=l2 ? l1-l2 : l2-l1)][1] = index++;
            }
        }
    }
    return lookup_table;
}

void c_destroy_bispectrum_lookup_table(c_blt table, size_t bandlimit)
{
    for (long l1 = 0; l1<=bandlimit; ++l1)
    {
        for (long l2=0; l2<=l1; ++l2)
        {
            for (long l = (l1>=l2 ? l1-l2 : l2-l1); l<=l1+l2 && l<=bandlimit; ++l)
            {
                free(table[l1][l2][l - (l1>=l2 ? l1-l2 : l2-l1)]);
            }
            free(table[l1][l2]);
        }
        free(table[l1]);
    }
    free(table);
}


r_bispectrum_lookout_table r_build_bispectrum_lookup_table(size_t bandlimit)
{
    size_t index = 0; 

    size_t ***lookup_table = malloc((bandlimit+1)*sizeof(size_t **));
    for (long l1 = 0; l1<=bandlimit; ++l1)
    {
        lookup_table[l1] = malloc((l1+1)*sizeof(size_t *));
        for (long l2=0; l2<=l1; ++l2)
        {
            lookup_table[l1][l2] = malloc((((l1+l2)>=bandlimit ? bandlimit : l1+l2)  - (l1>=l2 ? l1-l2 : l2-l1) + 1)*sizeof(size_t));
            for (long l = (l1>=l2 ? l1-l2 : l2-l1); l<=l1+l2 && l<=bandlimit; ++l)
            {
                lookup_table[l1][l2][l - (l1>=l2 ? l1-l2 : l2-l1)] = index;
                ++index;
            }
        }
    }
    return lookup_table;
}


void r_destroy_bispectrum_lookup_table(r_bispectrum_lookout_table table, size_t bandlimit)
{
    for (long l1 = 0; l1<=bandlimit; ++l1)
    {
        for (long l2=0; l2<=l1; ++l2)
        {
            free(table[l1][l2]);
        }
        free(table[l1]);
    }
    free(table);
}


/* Power spectrum fucntions */
double *c_allocate_power_spectrum(const size_t bandlimit)
{
    return malloc((bandlimit+1)*sizeof(double));
}


double *r_allocate_power_spectrum(const size_t bandlimit)
{
    return malloc((bandlimit+1)*sizeof(double));
}


void c_power_spectrum(c_shc * const shc, double *c_power_spectrum)
{
    double real_part; 
    double imag_part;

    for (long l = 0; l<=shc->bandlimit; ++l)
    {
        c_power_spectrum[l] = 0;
        for (long m = -l; m<=l; ++m)
        {
            real_part = c_get_shc(shc, REAL_PART, l, m);
            imag_part = c_get_shc(shc, IMAG_PART, l, m);
            c_power_spectrum[l] += real_part*real_part + imag_part*imag_part;
        }
    }
}


void r_power_spectrum(r_shc * const shc, double *r_pow_spec)
{
    double real_part;
    double imag_part;

    for (long l = 0; l<=shc->bandlimit; ++l)
    {
        real_part = r_get_shc(shc, REAL_PART, l, 0);
        r_pow_spec[l] = real_part*real_part;

        for (long m = 1; m<=l; ++m)
        {
            real_part = r_get_shc(shc, REAL_PART, l, m);
            imag_part = r_get_shc(shc, IMAG_PART, l, m);
            r_pow_spec[l] += 2.0*(real_part*real_part + imag_part*imag_part);
        }
    }
}


/* Bispectrum functions */

void c_bispectrum(c_shc * const shc, const c_blt lookup, cg_table * const table, double *c_bisp)
{
    // #pragma omp parallel for collapse(3)
    
    for (long l1 = 0; l1<=shc->bandlimit; ++l1)
    {
        for (long l2 = 0; l2<=l1; ++l2)
        {
            for (long l = l1-l2; l<=l1+l2 && l<=shc->bandlimit; ++l)
            {
                c_bisp[lookup[l1][l2][l-(l1-l2)][REAL_PART]]
                    = c_bispectral_invariant_real_part(shc, l1, l2, l, table);
                
                c_bisp[lookup[l1][l2][l-(l1-l2)][IMAG_PART]] 
                    = c_bispectral_invariant_imaginary_part(shc, l1, l2, l, table);
            }
        }
    }
}


double r_bispectral_invariant_real_part(r_shc * const shc, const long l1, const long l2, const long l, cg_table * const table)
{
    double invariant = 0;

    double U = 0;
    double M = 0;

    long lower_bound = 0;
    long upper_bound = 0;
    for (long m = -l; m<=l; ++m)
    {
        U = 0;
        M = 0;

        lower_bound = (-l1)>=(m-l2) ? (-l1) : (m-l2);
        upper_bound = l1>=(m+l2) ? (m+l2) : (l1);
        for (long m1 = lower_bound; m1<=upper_bound; ++m1)
        {
            U   +=  get_cg(table, l1, l2, l, m, m1)
                    *(
                        r_get_shc(shc, REAL_PART, l1, m1)*r_get_shc(shc, REAL_PART, l2, m-m1)
                            - r_get_shc(shc, IMAG_PART, l1, m1)*r_get_shc(shc, IMAG_PART, l2, m-m1)
                    );
                    
            M   +=  get_cg(table, l1, l2, l, m, m1)
                    *(
                        r_get_shc(shc, IMAG_PART, l1, m1)*r_get_shc(shc, REAL_PART, l2, m-m1)
                            + r_get_shc(shc, REAL_PART, l1, m1)*r_get_shc(shc, IMAG_PART, l2, m-m1)
                    );
        }

        invariant   +=  r_get_shc(shc, REAL_PART, l, m)*U 
                        + r_get_shc(shc, IMAG_PART, l, m)*M;
    }

    return invariant;
}


double r_bispectral_invariant_imaginary_part(r_shc * const shc, const long l1, const long l2, const long l, cg_table * const table)
{
    double invariant = 0;

    double U = 0;
    double M = 0;

    long lower_bound = 0;
    long upper_bound = 0;
    for (long m = -l; m<=l; ++m)
    {
        U = 0;
        M = 0;

        lower_bound = (-l1)>=(m-l2) ? (-l1) : (m-l2);
        upper_bound = l1>=(m+l2) ? (m+l2) : (l1);
        for (long m1 = lower_bound; m1<=upper_bound; ++m1)
        {
            U   += get_cg(table, l1, l2, l, m, m1)
                    *( 
                        r_get_shc(shc, REAL_PART, l1, m1)*r_get_shc(shc, REAL_PART, l2, m-m1)
                            - r_get_shc(shc, IMAG_PART, l1, m1)*r_get_shc(shc, IMAG_PART, l2, m-m1)
                    );

            M +=    get_cg(table, l1, l2, l, m, m1)
                    *(
                        r_get_shc(shc, IMAG_PART, l1, m1)*r_get_shc(shc, REAL_PART, l2, m-m1)
                        +   r_get_shc(shc, REAL_PART, l1, m1)*r_get_shc(shc, IMAG_PART, l2, m-m1)
                    );
        }

        invariant   += r_get_shc(shc, IMAG_PART, l, m)*U 
                        - r_get_shc(shc, REAL_PART, l, m)*M;
    }

    return invariant;
}



double c_bispectral_invariant_real_part(c_shc * const shc, const long l1, const long l2, const long l, const cg_table *table)
{
    double invariant = 0;

    double U = 0;
    double M = 0;

    long lower_bound = 0;
    long upper_bound = 0;
    for (long m = -l; m<=l; ++m)
    {
        U = 0;
        M = 0;

        lower_bound = (-l1)>=(m-l2) ? (-l1) : (m-l2);
        upper_bound = l1>=(m+l2) ? (m+l2) : (l1);
        for (long m1 = lower_bound; m1<=upper_bound; ++m1)
        {
            U   +=  get_cg(table, l1, l2, l, m, m1)
                    *(
                        c_get_shc(shc, REAL_PART, l1, m1)*c_get_shc(shc, REAL_PART, l2, m-m1)
                            - c_get_shc(shc, IMAG_PART, l1, m1)*c_get_shc(shc, IMAG_PART, l2, m-m1)
                    );
                    
            M   +=  get_cg(table, l1, l2, l, m, m1)
                    *(
                        c_get_shc(shc, IMAG_PART, l1, m1)*c_get_shc(shc, REAL_PART, l2, m-m1)
                            + c_get_shc(shc, REAL_PART, l1, m1)*c_get_shc(shc, IMAG_PART, l2, m-m1)
                    );
        }

        invariant   +=  c_get_shc(shc, REAL_PART, l, m)*U 
                        + c_get_shc(shc, IMAG_PART, l, m)*M;
    }

    return invariant;
}


double c_bispectral_invariant_imaginary_part(c_shc * const shc, const long l1, const long l2, const long l, const cg_table *table)
{
    double invariant = 0;

    double U = 0;
    double M = 0;

    long lower_bound = 0;
    long upper_bound = 0;
    for (long m = -l; m<=l; ++m)
    {
        U = 0;
        M = 0;

        lower_bound = (-l1)>=(m-l2) ? (-l1) : (m-l2);
        upper_bound = l1>=(m+l2) ? (m+l2) : (l1);
        for (long m1 = lower_bound; m1<=upper_bound; ++m1)
        {
            U   += get_cg(table, l1, l2, l, m, m1)
                    *( 
                        c_get_shc(shc, REAL_PART, l1, m1)*c_get_shc(shc, REAL_PART, l2, m-m1)
                            - c_get_shc(shc, IMAG_PART, l1, m1)*c_get_shc(shc, IMAG_PART, l2, m-m1)
                    );

            M +=    get_cg(table, l1, l2, l, m, m1)
                    *(
                        c_get_shc(shc, IMAG_PART, l1, m1)*c_get_shc(shc, REAL_PART, l2, m-m1)
                        +   c_get_shc(shc, REAL_PART, l1, m1)*c_get_shc(shc, IMAG_PART, l2, m-m1)
                    );
        }

        invariant   += c_get_shc(shc, IMAG_PART, l, m)*U 
                        - c_get_shc(shc, REAL_PART, l, m)*M;
    }

    return invariant;
}


void r_bispectrum(r_shc * const shc, const r_blt lookup, cg_table * const table, double *r_bisp)
{
    // #pragma omp parallel for collapse(3)
    
    for (long l1 = 0; l1<=shc->bandlimit; ++l1)
    {
        for (long l2 = 0; l2<=l1; ++l2)
        {
            for (long l = l1-l2; l<=l1+l2 && l<=shc->bandlimit; ++l)
            {
                size_t row_index = lookup[l1][l2][l-(l1-l2)];

                r_bisp[row_index] = r_bispectral_invariant_real_part(shc, l1, l2, l, table);
            }
        }
    }
}


void c_bispectrum_gradient(c_shc * const shc, const c_blt lookup, cg_table * const table, double *c_bisp_grad)
{
    // TODO
}

void r_bispectrum_gradient(r_shc * const shc, const r_blt lookup, cg_table * const table, double *r_bisp_grad)
{
    /* I assume the array r_bipsectrum_gradient is initialized with zeros.
    TODO: document the hell out of this, including an external file explaining the derivation. */
    for (long l1 = 0; l1<=shc->bandlimit; ++l1)
    {
        for (long l2 = 0; l2<=l1; ++l2)
        {
            for (long l = l1-l2; l<=l1+l2 && l<=shc->bandlimit; ++l)
            {
                double Uplus = 0;
                double Uminus = 0;
                double Mplus = 0;
                double Mminus = 0;

                // Working on the derivative of b_{l1,l2,l}
                size_t row_index = lookup[l1][l2][l-(l1-l2)]*(shc->bandlimit+1)*(shc->bandlimit+1);

                // *********************************
                // Derivative w.r.t f_{sigma,l,m}
                // *********************************
                long lower_bound = 0;
                long upper_bound = 0;
                
                for (long m = 1; m<=l; ++m) // for m>=1
                {
                    Uplus = 0;
                    Uminus = 0;
                    Mplus = 0;
                    Mminus = 0;

                    lower_bound = (-l1)>=(m-l2) ? (-l1) : (m-l2);
                    upper_bound = l1>=(m+l2) ? (m+l2) : (l1);
                    for (long m1 = lower_bound; m1<=upper_bound; ++m1)
                    {
                        Uplus   += get_cg(table, l1, l2, l, m, m1)
                                *( 
                                    r_get_shc(shc, REAL_PART, l1, m1)*r_get_shc(shc, REAL_PART, l2, m-m1)
                                    - r_get_shc(shc, IMAG_PART, l1, m1)*r_get_shc(shc, IMAG_PART, l2, m-m1)
                                );
                        Mplus   += get_cg(table, l1, l2, l, m, m1)
                                *(
                                    r_get_shc(shc, IMAG_PART, l1, m1)*r_get_shc(shc, REAL_PART, l2, m-m1)
                                    + r_get_shc(shc, REAL_PART, l1, m1)*r_get_shc(shc, IMAG_PART, l2, m-m1)
                                );
                    }

                    lower_bound = (-l1)>=(-m-l2) ? (-l1) : (-m-l2);
                    upper_bound = l1>=(-m+l2) ? (-m+l2) : (l1);
                    for (long m1 = lower_bound; m1<=upper_bound; ++m1)
                    {
                        Uminus  += get_cg(table, l1, l2, l, -m, m1)
                                *( 
                                    r_get_shc(shc, REAL_PART, l1, m1)*r_get_shc(shc, REAL_PART, l2, -m-m1)
                                    - r_get_shc(shc, IMAG_PART, l1, m1)*r_get_shc(shc, IMAG_PART, l2, -m-m1)
                                );
                        Mminus  += get_cg(table, l1, l2, l, -m, m1)
                                *(
                                    r_get_shc(shc, IMAG_PART, l1, m1)*r_get_shc(shc, REAL_PART, l2, -m-m1)
                                    + r_get_shc(shc, REAL_PART, l1, m1)*r_get_shc(shc, IMAG_PART, l2, -m-m1)
                                );
                    }

                    // The derivative w.r.t f_{r,l,m} (m>=1)
                    r_bisp_grad[row_index+l*l+2*m-1] += Uplus + ((m % 2)==0 ? 1 : -1)*Uminus;

                    // The derivative w.r.t f_{i,l,m} (m>=1)
                    r_bisp_grad[row_index+l*l+2*m] += Mplus - ((m % 2)==0 ? 1 : -1)*Mminus;
                }

                // The derivative w.r.t f_{r,l,0}
                lower_bound = (-l1)>=(-l2) ? (-l1) : (-l2);
                upper_bound = l1>=l2 ? (l2) : (l1);
                for (long m1 = lower_bound; m1<=upper_bound; ++m1)
                {
                    r_bisp_grad[row_index+l*l]  += get_cg(table, l1, l2, l, 0, m1)
                                                    *( 
                                                        r_get_shc(shc, REAL_PART, l1, m1)*r_get_shc(shc, REAL_PART, l2, -m1)
                                                        - r_get_shc(shc, IMAG_PART, l1, m1)*r_get_shc(shc, IMAG_PART, l2, -m1)
                                                    );
                }
                

                // *********************************
                // Derivative w.r.t f_{sigma,l1,m1}
                // *********************************
                lower_bound = 0;
                upper_bound = 0;
                
                for (long m1 = 1; m1<=l1; ++m1) // for m1>=1
                {
                    Uplus = 0;
                    Uminus = 0;
                    Mplus = 0;
                    Mminus = 0;

                    lower_bound = (-l)>=(m1-l2) ? (-l) : (m1-l2);
                    upper_bound = l>=(m1+l2) ? (m1+l2) : l;
                    for (long m = lower_bound; m<=upper_bound; ++m)
                    {
                        Uplus   += get_cg(table, l1, l2, l, m, m1)
                                *( 
                                    r_get_shc(shc, REAL_PART, l, m)*r_get_shc(shc, REAL_PART, l2, m-m1)
                                    + r_get_shc(shc, IMAG_PART, l, m)*r_get_shc(shc, IMAG_PART, l2, m-m1)
                                );

                        Mplus   += get_cg(table, l1, l2, l, m, m1)
                                *(
                                    r_get_shc(shc, IMAG_PART, l, m)*r_get_shc(shc, REAL_PART, l2, m-m1)
                                    - r_get_shc(shc, REAL_PART, l, m)*r_get_shc(shc, IMAG_PART, l2, m-m1)
                                );
                    }

                    lower_bound = (-l)>=(-m1-l2) ? (-l) : (-m1-l2);
                    upper_bound = l>=(-m1+l2) ? (-m1+l2) : l;
                    for (long m = lower_bound; m<=upper_bound; ++m)
                    {
                        Uminus  += get_cg(table, l1, l2, l, m, -m1)
                                *( 
                                    r_get_shc(shc, REAL_PART, l, m)*r_get_shc(shc, REAL_PART, l2, m+m1)
                                    + r_get_shc(shc, IMAG_PART, l, m)*r_get_shc(shc, IMAG_PART, l2, m+m1)
                                );

                        Mminus  += get_cg(table, l1, l2, l, m, -m1)
                                *(
                                    r_get_shc(shc, IMAG_PART, l, m)*r_get_shc(shc, REAL_PART, l2, m+m1)
                                    - r_get_shc(shc, REAL_PART, l, m)*r_get_shc(shc, IMAG_PART, l2, m+m1)
                                );
                    }

                    // The derivative w.r.t f_{r,l1,m1} (m1>=1)
                    r_bisp_grad[row_index+l1*l1+2*m1-1] += Uplus + ((m1 % 2)==0 ? 1 : -1)*Uminus;
                
                    // The derivative w.r.t f_{i,l1,m1} (m1>=1)
                    r_bisp_grad[row_index+l1*l1+2*m1] += Mplus - ((m1 % 2)==0 ? 1 : -1)*Mminus;
                }

                // The derivative w.r.t f_{r,l1,0}
                lower_bound = (-l)>=(-l2) ? (-l) : (-l2);
                upper_bound = l>=l2 ? l2 : l;
                for (long m = lower_bound; m<=upper_bound; ++m)
                {
                    r_bisp_grad[row_index+l1*l1]    += get_cg(table, l1, l2, l, m, 0)
                                                    *( 
                                                        r_get_shc(shc, REAL_PART, l, m)*r_get_shc(shc, REAL_PART, l2, m)
                                                        + r_get_shc(shc, IMAG_PART, l, m)*r_get_shc(shc, IMAG_PART, l2, m)
                                                    );
                }
                


                // *********************************
                // Derivative w.r.t f_{sigma,l2,m2}
                // *********************************
                lower_bound = 0;
                upper_bound = 0;
                
                for (long m2 = 1; m2<=l2; ++m2) // for m2>=1
                {
                    Uplus = 0;
                    Uminus = 0;
                    Mplus = 0;
                    Mminus = 0;

                    lower_bound = (-l)>=(m2-l1) ? (-l) : (m2-l1);
                    upper_bound = l>=(m2+l1) ? (m2+l1) : l;
                    for (long m = lower_bound; m<=upper_bound; ++m)
                    {
                        Uplus   += get_cg(table, l1, l2, l, m, m-m2)
                                *( 
                                    r_get_shc(shc, REAL_PART, l, m)*r_get_shc(shc, REAL_PART, l1, m-m2)
                                    + (shc, IMAG_PART, l, m)*r_get_shc(shc, IMAG_PART, l1, m-m2)
                                );

                        Mplus   += get_cg(table, l1, l2, l, m, m-m2)
                                *(
                                    r_get_shc(shc, IMAG_PART, l, m)*r_get_shc(shc, REAL_PART, l1, m-m2)
                                    - r_get_shc(shc, REAL_PART, l, m)*r_get_shc(shc, IMAG_PART, l1, m-m2)
                                );
                    }

                    lower_bound = (-l)>=(-m2-l1) ? (-l) : (-m2-l1);
                    upper_bound = l>=(-m2+l1) ? (-m2+l1) : l;
                    for (long m = lower_bound; m<=upper_bound; ++m)
                    {
                        Uminus  += get_cg(table, l1, l2, l, m, m+m2)
                                *( 
                                    r_get_shc(shc, REAL_PART, l, m)*r_get_shc(shc, REAL_PART, l1, m+m2)
                                    + r_get_shc(shc, IMAG_PART, l, m)*r_get_shc(shc, IMAG_PART, l1, m+m2)
                                );

                        Mminus  += get_cg(table, l1, l2, l, m, m+m2)
                                *(
                                    r_get_shc(shc, IMAG_PART, l, m)*r_get_shc(shc, REAL_PART, l1, m+m2)
                                    - r_get_shc(shc, REAL_PART, l, m)*r_get_shc(shc, IMAG_PART, l1, m+m2)
                                );
                    }

                    // The derivative w.r.t f_{r,l2,m2} (m2>=1)
                    r_bisp_grad[row_index+l2*l2+2*m2-1] += Uplus + ((m2 % 2)==0 ? 1 : -1)*Uminus;
                
                    // The derivative w.r.t f_{i,l2,m2} (m2>=1)
                    r_bisp_grad[row_index+l2*l2+2*m2] += Mplus - ((m2 % 2)==0 ? 1 : -1)*Mminus;
                }

                // The derivative w.r.t f_{r,l2,0}
                lower_bound = (-l)>=(-l1) ? (-l) : (-l1);
                upper_bound = l>=l1 ? l1 : l;
                for (long m = lower_bound; m<=upper_bound; ++m)
                {
                    r_bisp_grad[row_index+l2*l2]    += get_cg(table, l1, l2, l, m, m)
                                                    *( 
                                                        r_get_shc(shc, REAL_PART, l, m)*r_get_shc(shc, REAL_PART, l1, m)
                                                        + r_get_shc(shc, IMAG_PART, l, m)*r_get_shc(shc, IMAG_PART, l1, m)
                                                    );
                }
            }
        }
    }
}


/* Printing functions */
void r_print_power_spectrum(double *r_pow_spec, size_t bandlimit)
{
    printf("Bandlimit = %d.\n", bandlimit);
    printf("=========================================================================\n");
    printf("\tOrder\t\tValue\n");
    printf("=========================================================================\n");
    for (long l = 0; l<=bandlimit; ++l)
    {
        printf("\t%3d\t\t% 14.14f\n", l, r_pow_spec[l]);
    }
}


void c_print_power_spectrum(double *c_pow_spec, size_t bandlimit)
{
    r_print_power_spectrum(c_pow_spec, bandlimit);
}


void r_print_all_bispectral_invariants(r_shc * const shc, cg_table * const table)
{
    printf("Bandlimit = %d.\n", shc->bandlimit);
    printf("=========================================================================\n");
    printf("\t( l1, l2, l)\t\tReal part\t\tImag part\n");
    printf("=========================================================================\n");
    for (long l1 = 0; l1<=shc->bandlimit; ++l1)
    {
        for (long l2 = 0; l2<=l1; ++l2)
        {
            for (long l = l1-l2; l<=shc->bandlimit && l<=l1+l2; ++l)
            {
                printf("\t(% 3d,% 3d,% 3d)\t\t% 14.14f\t% 14.14f\n", 
                        l1, l2, l, 
                        r_bispectral_invariant_real_part(shc, l1, l2, l, table), 
                        r_bispectral_invariant_imaginary_part(shc, l1, l2, l, table));
                // printf("\t(% 3d,% 3d,% 3d)\t\t%14.14f\n", 
                        // l1, l2, l, r_bispectral_invariant_imaginary_part(shc, l1, l2, l, table));
            }
        }
    }
}


void c_print_all_bispectral_invariants()
{
    // TODO
}
