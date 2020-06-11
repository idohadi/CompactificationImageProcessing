// TODO: write docs in doxygen

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "spherical_harmonics.h"

int main(int argc, char *argv[])
{
    sfmt_t sfmt;
    sfmt_init_gen_rand(&sfmt, time(NULL));

    if (argc==3)
    {
        size_t bandlimit;
        sscanf(argv[2], "%zu", &bandlimit);

        if (argv[1][0]=='c')
        {
            c_shc shc = c_init_shc2(bandlimit);
            c_random_normalized_shc(&sfmt, &shc);
            c_print_shc(&shc);

            c_deallocate_coefficients(shc.coefficients);
            printf("Done with everything, including deallocating.");
            return 0;
        }
        if (argv[1][0]=='r')
        {
            r_shc shc = r_init_shc2(bandlimit);
            r_random_normalized_shc(&sfmt, &shc);
            r_print_shc(&shc);

            r_deallocate_coefficients(shc.coefficients);
            printf("Done with everything, including deallocating.");
            return 0;
        }
        
        printf("Can't understand arguments.");
        return 0;
    }
    else
    {
        printf("Wrong number of arguments.");
    }
}
