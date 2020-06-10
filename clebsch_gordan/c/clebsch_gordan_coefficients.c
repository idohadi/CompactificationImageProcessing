// TODO: write docs in doxygen


/* Function declarations. Defined below */
#include "clebsch_gordan_coefficients.h"


/* Clebsch-Gordan management functions */
long cg_lower_bound(const long l1, const long l2, const long m)
{
    /* Returns the lowest m-index of Clebsch-Gordan coefficients of order (l1, l2, l, m)  */

    return (m-l2<=-l1) ? -l1 : m-l2;
}


long cg_upper_bound(const long l1, const long l2, const long m)
{
    /* Returns the highest m-index of Clebsch-Gordan coefficients of order (l1, l2, l, m)  */

    return (m+l2<=l1) ? m+l2 : l1;
}


double backward_substitution_central_diagonal(const long l1, const long l2, const long l, const long m1, const long m)
{
    return l1*(l1+1) + l2*(l2+1) + 2*m1*(m-m1) - l*(l+1);
}

double backward_substitution_secondary_diagonal(const long l1, const long l2, const long m1, const long m)
{
    return sqrt(l1*(l1+1) - m1*(m1+1)) * sqrt(l2*(l2+1) - (m-m1)*(m-m1-1));
}


void calculate_clebsch_gordan(const long l1, const long l2, const long l, const long m, double * restrict cg)
{
    /*  I assume cg is an array of size (upper_bound - lower_bound + 1). */
    long lower_bound = clebsch_gordan_lower_bound(l1, l2, m);
    long upper_bound = clebsch_gordan_upper_bound(l1, l2, m);

    cg[upper_bound - lower_bound] = 1.0/(upper_bound - lower_bound + 1);

    if (upper_bound > lower_bound)
    {   
        double norm = cg[upper_bound - lower_bound]*cg[upper_bound - lower_bound];

        // Using the recursion to calculate the Clebsch-Gordan coefficients
        cg[upper_bound - lower_bound - 1] 
            = -cg[upper_bound - lower_bound]*backward_substitution_central_diagonal(l1, l2, l, upper_bound, m)/backward_substitution_secondary_diagonal(l1, l2, upper_bound-1, m);
        norm += cg[upper_bound - lower_bound - 1]*cg[upper_bound - lower_bound - 1];

        for (long i = upper_bound - lower_bound - 2; i>=0; --i)
        {
            cg[i] = - (cg[i+1]*backward_substitution_central_diagonal(l1, l2, l, lower_bound+i+1, m)
                     + cg[i+2]*backward_substitution_secondary_diagonal(l1, l2, lower_bound+i+1, m))
                              /backward_substitution_secondary_diagonal(l1, l2, lower_bound+i, m);
            norm += cg[i]*cg[i];
        }
        
        // Normalize the solution so that the norm of the Clebsch-Gordan coefficients is 1
        norm = sqrt(norm);
        for (long i = 0; i<upper_bound-lower_bound+1; ++i)
        {
            cg[i] /= norm;
        }
    }
}


void calculate_clebsch_gordan_coefficients(const size_t bandlimit, cg_table clebsch_gordan_table)
{
    /* Documentation of memory structure is in allocate_memeory_for_clebsch_gordan_coefficients */

    #pragma omp parallel for collapse(2)
    for (long l1 = 0; l1<=bandlimit; ++l1)
    {
        for (long l2 = 0; l2<=bandlimit; ++l2)
        {
            for (long l = absolute_value(l1-l2); l<=minimum(l1+l2, bandlimit); ++l)
            {
                for (long m = -l; m<=l; ++m)
                {
                    calculate_clebsch_gordan(l1, l2, l, m, clebsch_gordan_table[l1][l2][l-absolute_value(l1-l2)][m+l]);
                }
            }
        }
    }
}


cg_table allocate_memeory_for_clebsch_gordan_coefficients(const size_t bandlimit)
{
    /* 
        Returns an array CGs with dimensions such that it is possible to assign:
            <l1 m1 l2 (m-m1) | l m> = CGs[l1][l2][l-abs(l1-l2)][m+l][m1-min_m1]
        where
            0<=l1<=bandlimit
            0<=l2<=bandlimit
            abs(l1-l2)<=l<=min(l1+l2, bandlimit)
            -l<=m<=l
            min_m1<=m1<=max_m1
                min_m1:=max(-l1, m-l2)
                max_m1:=min(l1,  m+l2)  */

    cg_table clebsch_gordan_table = malloc((bandlimit+1)*sizeof(double ****));

    for (long l1 = 0; l1<=bandlimit; ++l1)
    {
        clebsch_gordan_table[l1] = malloc((bandlimit+1)*sizeof(double ***));
        
        for (long l2 = 0; l2<=bandlimit; ++l2)
        {
            clebsch_gordan_table[l1][l2] = malloc((minimum(l1+l2, bandlimit) - absolute_value(l1-l2) + 1)*sizeof(double **));

            for (long l = absolute_value(l1-l2); l<=minimum(l1+l2, bandlimit); ++l)
            {
                clebsch_gordan_table[l1][l2][l- absolute_value(l1-l2)] = malloc((2*l+1)*sizeof(double *));

                for (long m = -l; m<=l; ++m)
                {
                    clebsch_gordan_table[l1][l2][l- absolute_value(l1-l2)][m+l] 
                        = malloc((clebsch_gordan_upper_bound(l1, l2, m) - clebsch_gordan_lower_bound(l1, l2, m) + 1)*sizeof(double));
                }
            }
        }
    }

    return clebsch_gordan_table;
}

void free_memory_for_clebsch_gordan_coefficients(const size_t bandlimit, cg_table clebsch_gordan_table)
{
    /* Documentation of memory structure is in allocate_memeory_for_clebsch_gordan_coefficients */

    for (long l1 = 0; l1<=bandlimit; ++l1)
    {
        for (long l2 = 0; l2<=bandlimit; ++l2)
        {
            for (long l = absolute_value(l1-l2); l<=minimum(l1+l2, bandlimit); ++l)
            {
                for (long m = -l; m<=l; ++m)
                {
                    free(clebsch_gordan_table[l1][l2][l- absolute_value(l1-l2)][m+l]);
                }

                free(clebsch_gordan_table[l1][l2][l- absolute_value(l1-l2)]);
            }

            free(clebsch_gordan_table[l1][l2]);
        }

        free(clebsch_gordan_table[l1]);
    }

    free(clebsch_gordan_table);
}


bool validate_clebsch_gordan_table_order(const size_t bandlimit, const long l1, const long l2, const long l, const long m, const long m1)
{
    /* Documentation of memory structure is in allocate_memeory_for_clebsch_gordan_coefficients */
    
    if (l1 > bandlimit || l1 < 0)
    {
        return false;
    }

    if (l2 > bandlimit || l2 < 0)
    {
        return false;
    }

    if (l < absolute_value(l1-l2) || l > minimum(l1+l2, bandlimit))
    {
        return false;
    }

    if (m < -l || m > l)
    {
        return false;
    }

    if (m < clebsch_gordan_lower_bound(l1, l2, m) || m > clebsch_gordan_upper_bound(l1, l2, m))
    {
        return false;
    }

    return true;
}


double get_cg(const cg_table table, const long l1, const long l2, const long l, const long m, const long m1)
{
    /* Documentation of memory structure is in allocate_memeory_for_clebsch_gordan_coefficients */
    
    return table[l1][l2][l-absolute_value(l1-l2)][m+l][m1-clebsch_gordan_lower_bound(l1, l2, m)];
}


void set_cg(cg_table clebsch_gordan_table, const long l1, const long l2, const long l, const long m, const long m1, const double coefficient)
{
    /* Documentation of memory structure is in allocate_memeory_for_clebsch_gordan_coefficients */
    
    clebsch_gordan_table[l1][l2][l-absolute_value(l1-l2)][m+l][m1-clebsch_gordan_lower_bound(l1, l2, m)] 
        = coefficient;
}
