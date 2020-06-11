// TODO: write docs in doxygen
/* 
USE:
        spherical_spectra_tester flag bandlimit flag2
    where 
        flag = c for complex-valued spherical function, r for real-valued spherical function
        bandlimit is non-negative integer
        flag2 = p for power spectrum, b for bispectrum
*/

#include <time.h>
#include "spherical_spectra.h"

int main(int argc, char *argv[])
{
    sfmt_t sfmt;
    sfmt_init_gen_rand(&sfmt, time(NULL));

    if (argc==4)
    {
        size_t bandlimit;
        sscanf(argv[2], "%zu", &bandlimit);

        if (argv[1][0]=='c')
        {
            c_shc shc = c_init_shc2(bandlimit);
            c_random_normalized_shc(&sfmt, &shc);
            c_print_shc(&shc);

            if (argv[3][0]=='p')
            {
                double *powspec = c_allocate_power_spectrum(bandlimit);
                c_power_spectrum(&shc, powspec);

                c_print_power_spectrum(powspec, bandlimit);
            }
            
            if (argv[3][0]=='b')
            {
                printf("TODO: Calculate and print the bipsectrum for complex-valued function.\n");
            }

            c_deallocate_coefficients(shc.coefficients);
            printf("Done with everything, including deallocating.");

            return 0;
        }
        if (argv[1][0]=='r')
        {
            r_shc shc = r_init_shc2(bandlimit);
            r_random_normalized_shc(&sfmt, &shc);
            r_print_shc(&shc);

            if (argv[3][0]=='p')
            {
                double *powspec = r_allocate_power_spectrum(bandlimit);
                r_power_spectrum(&shc, powspec);

                r_print_power_spectrum(powspec, bandlimit);
            }
            
            if (argv[3][0]=='b')
            {
                cg_table table = allocate_cg_table(bandlimit);
                calculate_cg_table(bandlimit, &table);

                r_print_all_bispectral_invariants(&shc, &table);
                destroy_cg_table(&table);
            }
            
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
