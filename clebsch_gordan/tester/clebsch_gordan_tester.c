// TODO: write docs in doxygen

#include <stdio.h>
#include "clebsch_gordan_coefficients.h"

int main()
{
    printf("Enter l1, l2, l, m:\n");
    long l1, l2, l, m;
    scanf("%d, $d, $d, $d", l1, l2, l, m);
    print_cg_vector(l1, l2, l, m);
    
    return 0;
}
