#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include "clebsch_gordan_coefficients.h"
#include "spherical_spectra.h"


int main()
{
    size_t bandlimit = 3;
    size_t ***lookup_table2 = r_build_bispectrum_lookup_table(bandlimit);
    printf("Bandlimit = %d.\n", bandlimit);
    for (long int l1=0; l1<=bandlimit; ++l1)
    {
        for (long int l2=0; l2<=l1; ++l2)
        {
            for (long int l=absolute_value(l1-l2); l<=minimum(l1+l2, bandlimit); ++l)
            {
                printf("Lookup table value for (l1,l2,l)=(%d, %d, %d) is \t %d\n", 
                        l1, l2, l, lookup_table2[l1][l2][l-absolute_value(l1-l2)]+1);
            }
        }
    }

    // cb_table cb_coeffs = allocate_memeory_for_clebsch_gordan_coefficients(bandlimit);

    // double *grad = malloc(((bandlimit+1)*(bandlimit+1)*lookup_table2[bandlimit][bandlimit][bandlimit]+1)*sizeof(double));

    // double shc[81] = 
    // {-0.3448691318495896,0.3425287409034796,-0.1227100348260888,0.6670011911779505,0.5377085048592298,-0.6654929090105566,0.7239609574041443,0.9797443072630079,0.02884691301140885,0.7685620462539107,0.176052110616995,-0.6904953026879106,-0.6002743542850959,-0.1860903257221864,0.4974114364313829,0.6511676315723112,0.5799260598890617,-0.3629515092020157,0.06812825474145279,-0.8200986424588379,-0.7765885116135931,-0.7274149021234027,0.3573046096003765,-0.009645961820678783,-0.6205791879648397,-0.009988350019558334,-0.7047835560466227,-0.8900517061876236,0.7014253485780149,0.1211190547097696,0.8592177335133264,0.3933344011104551,0.1655819303516801,0.6307944229548423,0.7580278091943555,0.9778232321591782,-0.9989552492861105,0.7308771820260491,0.2251329389679975,0.9799004114176617,0.05536013867688472,-0.04095322957956249,0.6026952110439046,-0.5443141285879163,-0.003811417607220591,0.8017049770640097,0.1493224382603753,0.6903563701080733,0.4772805839908036,0.1719740716529516,-0.5065309480280502,0.332832434638936,-0.8330343727947547,0.2519195703431665,0.3218891158946846,0.4595037106344422,0.7815042326506445,0.9646064457672128,0.5380581706717924,0.1628929757507955,0.856626124628376,0.1601807315168833,-0.9660341233254774,-0.7582808578028837,0.7254214373993393,-0.03140697757579503,0.6897113491525264,-0.5811898319581306,0.10458268307755,0.2597667701288424,-0.9360179684748662,0.229426838234281,-0.2751770754538947,-0.9009348419158776,-0.0208600216453565,-0.6149792078758505,-0.7538325049081096,-0.5890116581846403,-0.7069701787702203,-0.6218556510547724,-0.9146951781777133};

    // r_calculate_bispectrum_gradient(shc, bandlimit, lookup_table2, cb_coeffs, grad);

    // for (long int l1=0; l1<=bandlimit; ++l1)
    // {
    //     for (long int l2=0; l2<=l1; ++l2)
    //     {
    //         for (long int l=absolute_value(l1-l2); l<=minimum(l1+l2, bandlimit); ++l)
    //         {
    //             for (long int l_der = 0; l_der<=bandlimit; ++l_der)
    //             {
    //                 for (long int m_der = 0; m_der<=l; ++m_der)
    //                 {
    //                     size_t row_index = lookup_table2[l1][l2][l-(l1-l2)]*(bandlimit+1)*(bandlimit+1);

    //                     if (m_der==0)
    //                     {
    //                         printf("gradient of b_{l1,l2,l}=b_{%d, %d, %d} w.r.t f_{r,l',m'}=f_{r,%d,%d} is \t %.10e\n", 
    //                             l1, l2, l, l_der, m_der, grad[row_index+l_der*l_der]);
    //                     }
    //                     else
    //                     {
    //                         printf("gradient of b_{l1,l2,l}=b_{%d, %d, %d} w.r.t f_{r,l',m'}=f_{r,%d,%d} is \t %.10e\n", 
    //                             l1, l2, l, l_der, m_der, grad[row_index+l_der*l_der+2*m_der-1]);
    //                         printf("gradient of b_{l1,l2,l}=b_{%d, %d, %d} w.r.t f_{i,l',m'}=f_{r,%d,%d} is \t %.10e\n", 
    //                             l1, l2, l, l_der, m_der, grad[row_index+l_der*l_der+2*m_der]);
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }
    // printf("\nCB Coeffs for bandlimit %d.\n", bandlimit);
    // cb_table cb_coeffs = allocate_memeory_for_clebsch_gordan_coefficients(bandlimit);
    // calculate_clebsch_gordan_coefficients(bandlimit, cb_coeffs);
    // double *c2; 

    // long int lb = 0;
    // long int ub = 0;

    // for (long int l1=0; l1<=bandlimit; ++l1)
    // {
    //     for (long int l2=0; l2<=l1; ++l2)
    //     {
    //         for (long int l=absolute_value(l1-l2); l<=minimum(l1+l2, bandlimit); ++l)
    //         {
    //             for (long int m = -l; m<=l; ++m)
    //             {
    //                 lb = clebsch_gordan_lower_bound(l1, l2, m);
    //                 ub = clebsch_gordan_upper_bound(l1, l2, m);

    //                 double *c1 = malloc((ub - lb + 1)*sizeof(double));
                    
    //                 for (long int i = 0; i<=(ub-lb+1); ++i)
    //                 {
    //                     c1[i] = 0;
    //                 }
                    
    //                 calculate_clebsch_gordan(l1, l2, l, m, c1);
                    
    //                 double norm = 0;
    //                 for (long int i = 0; i<=(ub-lb); ++i)
    //                 // for (long int i = lb; i<=ub; ++i)
    //                 {
    //                     norm += (c1[i] - get_clebsch_gordan_coefficient(cb_coeffs, l1, l2, l, m, lb+i))*(c1[i] - get_clebsch_gordan_coefficient(cb_coeffs, l1, l2, l, m, lb+i));
    //                     // norm += (get_clebsch_gordan_coefficient(cb_coeffs, l1, l2, l, m, i))*(get_clebsch_gordan_coefficient(cb_coeffs, l1, l2, l, m, i));
    //                 }
    //                 norm = sqrt(norm);

    //                 printf("Norm of (l1,l2,l,m)=(%d, %d, %d, %d) is \t %f\n", 
    //                     l1, l2, l, m, norm);
    //                 free(c1);
    //             }
                
    //         }
    //     }
    // }
}