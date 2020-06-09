#include <stdio.h>
#include <stdint.h>
#include "clebsch_gordan_coefficients.h"
#include "spherical_spectra.h"

int main()
{
    size_t bandlimit = 3;
    double f[16] = {   
          -0.804919190001181,
  -0.443003562265903,
   0.093763038409968,
   0.915013670868595,
   0.929777070398553,
  -0.684773836644903,
   0.941185563521231,
   0.914333896485891,
  -0.029248702554318,
   0.600560937777600,
  -0.716227322745569,
  -0.156477434747450,
   0.831471050378134,
   0.584414659119109,
   0.918984852785806,
   0.311481398313174,
   };


   double f2[16] ={   
       -0.804919190001181,
   0.903354705201618,
  -0.412793169215987,
  -0.604768171353372,
  -0.571675419393334,
   0.211309433977652,
   1.480016295128577,
  -0.169633600520792,
  -0.443360205463689,
   0.705889542568727,
   0.150377932268827,
  -0.874236336027483,
  -1.089626593808205,
  -0.354571388357020,
  -0.172275414753054,
  -0.559922663568557,
   }; // rotated version of the previous one

    cb_table cb_coeffs = allocate_memeory_for_clebsch_gordan_coefficients(bandlimit);
    calculate_clebsch_gordan_coefficients(bandlimit, cb_coeffs);

    for (long int l1 = 0; l1<=bandlimit; ++l1)
    {
        for (long int l2 = 0; l2<=l1; ++l2)
        {
            for (long int l = l1-l2; l<=l1+l2 && l<=bandlimit; ++l)
            {
                printf("Bispectral invariant (l1,l2,l)=(%d,%d,%d), real part = \t %f\n", 
                        l1, l2, l, 
                        r_calculate_bispectral_invariant_real_part(f, bandlimit, l1, l2, l, cb_coeffs));
                printf("Bispectral invariant (l1,l2,l)=(%d,%d,%d), imaginary part = \t %f\n", 
                        l1, l2, l, 
                        r_calculate_bispectral_invariant_imaginary_part(f, bandlimit, l1, l2, l, cb_coeffs));
            }
        }
    }
    // for (long int l = 0; l<=bandlimit; ++l)
    // {
    //     for (long int m = -l; m<=l; ++m)
    //     {
    //         printf("f_{r,%d,%d} = %f\t", l, m, r_get_spherical_harmonics(f, 'r', l, m));
    //         printf("f_{i,%d,%d} = %f\n", l, m, r_get_spherical_harmonics(f, 'i', l, m));
    //     }
    // }

    // size_t ***lookup_table = r_build_bispectrum_lookup_table(bandlimit);
    // cb_table cb_coeffs = allocate_memeory_for_clebsch_gordan_coefficients(bandlimit);
    // calculate_clebsch_gordan_coefficients(bandlimit, cb_coeffs);

    // printf("Finished perliminaries. Calculating bispectrum.\n");
    
    // printf("Largest index of bispectrum = %d\n", lookup_table[bandlimit][bandlimit][bandlimit]);

    // for (long int l1=0; l1<=bandlimit; ++l1)
    // {
    //     for (long int l2=0; l2<=l1; ++l2)
    //     {
    //         for (long int l=absolute_value(l1-l2); l<=minimum(l1+l2, bandlimit); ++l)
    //         {
    //             printf("Lookup table value for (l1,l2,l)=(%d, %d, %d) is \t %d\n", 
    //                     l1, l2, l, lookup_table[l1][l2][l-absolute_value(l1-l2)]);
    //         }
    //     }
    // }

    // printf("\n");

    // size_t bisp_size = lookup_table[bandlimit][bandlimit][bandlimit]+1;
    // double *bisp1 = malloc(bisp_size*sizeof(double));
    // double *bisp2 = malloc(bisp_size*sizeof(double));
    // for (int i = 0; i<bisp_size; ++i)
    // {
    //     bisp1[i] = 0.0;
    //     bisp2[i] = 0.0;
    // }

    // printf("Calculating bispectrum.\n");
    // r_calculate_bispectrum(f, bandlimit, lookup_table, cb_coeffs, bisp1);
    // r_calculate_bispectrum(f2, bandlimit, lookup_table, cb_coeffs, bisp2);
    // printf("Finished calculating bispectrum. Here it is:\n");
    // for (long int l1=0; l1<=bandlimit; ++l1)
    // {
    //     for (long int l2=0; l2<=l1; ++l2)
    //     {
    //         for (long int l=absolute_value(l1-l2); l<=minimum(l1+l2, bandlimit); ++l)
    //         {
    //             printf("(l1,l2,l)=(%d, %d, %d), ind = %d.\t  bisp1 =  %f\tbisp2 = %f\n", 
    //                     l1, l2, l, lookup_table[l1][l2][l-absolute_value(l1-l2)], 
    //                     bisp1[lookup_table[l1][l2][l-absolute_value(l1-l2)]], bisp2[lookup_table[l1][l2][l-absolute_value(l1-l2)]]);
    //         }
    //     }
    // }
    
    // bandlimit = 6;
    // size_t ***lookup_table2 = r_build_bispectrum_lookup_table(bandlimit);
    // printf("Bandlimit = %d.\n", bandlimit);
    // for (long int l1=0; l1<=bandlimit; ++l1)
    // {
    //     for (long int l2=0; l2<=l1; ++l2)
    //     {
    //         for (long int l=absolute_value(l1-l2); l<=minimum(l1+l2, bandlimit); ++l)
    //         {
    //             printf("Lookup table value for (l1,l2,l)=(%d, %d, %d) is \t %d\n", 
    //                     l1, l2, l, lookup_table2[l1][l2][l-absolute_value(l1-l2)]+1);
    //         }
    //     }
    // }
}