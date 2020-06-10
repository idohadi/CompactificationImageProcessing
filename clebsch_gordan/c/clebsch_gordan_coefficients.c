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
    long lower_bound = cg_lower_bound(l1, l2, m);
    long upper_bound = cg_upper_bound(l1, l2, m);

    cg[upper_bound - lower_bound] = 1.0/(upper_bound - lower_bound + 1);

    if (upper_bound > lower_bound)
    {   
        double norm = cg[upper_bound - lower_bound]*cg[upper_bound - lower_bound];

        // Using the recursion to calculate the Clebsch-Gordan coefficients
        cg[upper_bound - lower_bound - 1] 
            = -cg[upper_bound - lower_bound]*backward_substitution_central_diagonal(l1, l2, l, upper_bound, m)
                /backward_substitution_secondary_diagonal(l1, l2, upper_bound-1, m);
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
    /* Documentation of memory structure is in allocate_cg_table */

    #pragma omp parallel for collapse(2)
    for (long l1 = 0; l1<=bandlimit; ++l1)
    {
        for (long l2 = 0; l2<=bandlimit; ++l2)
        {
            for (long l = (l1>=l2 ? l1-l2 : l2-l1); l<=l1+l2 && l<=bandlimit; ++l)
            {
                for (long m = -l; m<=l; ++m)
                {
                    calculate_clebsch_gordan(l1, l2, l, m, clebsch_gordan_table[l1][l2][l-(l1>=l2 ? l1-l2 : l2-l1)][m+l]);
                }
            }
        }
    }
}


cg_table allocate_cg_table(const size_t bandlimit)
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

    cg_table table = malloc((bandlimit+1)*sizeof(double ****));

    long min_val; 
    long abs_val;

    for (long l1 = 0; l1<=bandlimit; ++l1)
    {
        table[l1] = malloc((bandlimit+1)*sizeof(double ***));
        
        for (long l2 = 0; l2<=bandlimit; ++l2)
        {
            min_val = l1+l2>=bandlimit ? bandlimit : l1+l2;
            abs_val = l1>=l2 ? l1-l2 : l2-l1;

            table[l1][l2] = malloc((min_val - abs_val + 1)*sizeof(double **));

            for (long l = abs_val; l<=min_val; ++l)
            {
                table[l1][l2][l - abs_val] = malloc((2*l+1)*sizeof(double *));

                for (long m = -l; m<=l; ++m)
                {
                    table[l1][l2][l - abs_val][m+l] 
                        = malloc((cg_upper_bound(l1, l2, m) - cg_lower_bound(l1, l2, m) + 1)*sizeof(double));
                }
            }
        }
    }

    return table;
}

void destroy_cg_table(const size_t bandlimit, cg_table table)
{
    /* Documentation of memory structure is in allocate_cg_table */

    long min_val; 
    long abs_val;

    for (long l1 = 0; l1<=bandlimit; ++l1)
    {
        for (long l2 = 0; l2<=bandlimit; ++l2)
        {
            min_val = l1+l2>=bandlimit ? bandlimit : l1+l2;
            abs_val = l1>=l2 ? l1-l2 : l2-l1;

            for (long l = abs_val; l<=min_val; ++l)
            {
                for (long m = -l; m<=l; ++m)
                {
                    free(table[l1][l2][l - abs_val][m+l]);
                }

                free(table[l1][l2][l- abs_val]);
            }

            free(table[l1][l2]);
        }

        free(table[l1]);
    }

    free(table);
}


double get_cg(const cg_table table, const long l1, const long l2, const long l, const long m, const long m1)
{
    /* Documentation of memory structure is in allocate_cg_table */
    
    return table[l1][l2][l-(l1>=l2 ? l1-l2 : l2-l1)][m+l][m1-cg_lower_bound(l1, l2, m)];
}


void set_cg(cg_table clebsch_gordan_table, const long l1, const long l2, const long l, const long m, const long m1, const double coefficient)
{
    /* Documentation of memory structure is in allocate_cg_table */
    
    clebsch_gordan_table[l1][l2][l-(l1>=l2 ? l1-l2 : l2-l1)][m+l][m1-cg_lower_bound(l1, l2, m)] 
        = coefficient;
}
