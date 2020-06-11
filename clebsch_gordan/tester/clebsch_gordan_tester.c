// TODO: write docs in doxygen

#include <stdio.h>
#include <stdlib.h>
#include "clebsch_gordan_coefficients.h"

int main(int argc, char *argv[])
{
    if (argc==5)
    {
        // Handling input
        long l1, l2, l, m;
        sscanf(argv[1], "%ld", &l1);
        sscanf(argv[2], "%ld", &l2);
        sscanf(argv[3], "%ld", &l);
        sscanf(argv[4], "%ld", &m);
        
        // Print output
        printf("\nInput: (l1,l2,l,m) = (%-3ld,%-3ld,%-3ld,%-3ld)\n", l1, l2, l, m);
        printf("Lower bound = %-3d. Upper bound = %-3d.\n\n", cg_lower_bound(l1, l2, m), cg_upper_bound(l1, l2, m));
        print_cg_vector(l1, l2, l, m);
    }
    else if (argc==2)
    {
        // Handle input
        size_t bandlimit;
        sscanf(argv[1], "%zu", &bandlimit);
        
        // Calculate Clebsch-Gordan table
        cg_table table = allocate_cg_table(bandlimit);
        calculate_cg_table(bandlimit, &table);

        // Print the table
        print_cg(&table);
    }
    else
    {
        printf("Wrong number of arguments.\n");
        return 0;
    }

    return 0;
}
