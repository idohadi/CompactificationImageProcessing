/** 
 * This is a tester. 
 * 
 * How to compile:
 *      TODO */

#include <stdio.h>
#include "clebsch_gordan_coefficients.h"

int main()
{
    size_t bandlimit = 3;
    
    printf("Allocating memory for bandlimit = %d\n", bandlimit);
    cb_table table = allocate_memeory_for_clebsch_gordan_coefficients(bandlimit);

    printf("Disallocating memory for bandlimit = %d\n\n", bandlimit);
    free_memory_for_clebsch_gordan_coefficients(bandlimit, table);
    printf("Disallocation successful.");

    bandlimit = 2;
    printf("Allocating memeory for bandlimit = %d\n", bandlimit);
    cb_table table2 = allocate_memeory_for_clebsch_gordan_coefficients(bandlimit);

    printf("Calculating coefficients for bandlimit = %d\n", bandlimit);
    calculate_clebsch_gordan_coefficients(bandlimit, table2);

    printf("Printing coefficients for bandlimit = %d\n", bandlimit);
        for (long int l1 = 0; l1<=bandlimit; ++l1)
    {
        for (long int l2 = 0; l2<=bandlimit; ++l2)
        {
            for (long int l = absolute_value(l1-l2); l<=minimum(l1+l2, bandlimit); ++l)
            {
                for (long int m = -l; m<=l; ++m)
                {
                    printf("l1=%d, l2=%d, l=%d, m=%d\n", l1, l2, l, m);
                    for (long int m1 = clebsch_gordan_lower_bound(l1, l2, m); m1<=clebsch_gordan_upper_bound(l1, l2, m); ++m1)
                    {
                        printf("Coeff %d = %f\n", m1, get_clebsch_gordan_coefficient(table2, l1, l2, l, m, m1));
                    }
                    printf("\n");
                }
            }
        }
    }
    
    // bandlimit = 50;
    // printf("Allocating memeory for bandlimit = %d\n", bandlimit);
    // cb_table table3 = allocate_memeory_for_clebsch_gordan_coefficients(bandlimit);
    // printf("Calculating coefficients for bandlimit = %d\n", bandlimit);
    // calculate_clebsch_gordan_coefficients(bandlimit, table3);
    // printf("Done calculating coefficients for bandlimit = %d\n", bandlimit);
}