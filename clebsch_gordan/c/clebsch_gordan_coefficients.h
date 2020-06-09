// TODO: make sure it's correct

#include <stdint.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

typedef double***** cg_table;

long int clebsch_gordan_lower_bound(const long int l1, const long int l2, const long int m);

long int clebsch_gordan_upper_bound(const long int l1, const long int l2, const long int m);

void calculate_clebsch_gordan(const long int l1, const long int l2, const long int l, const long int m, double *cb);

void calculate_clebsch_gordan_coefficients(const size_t bandlimit, cg_table clebsch_gordan_table);

cg_table allocate_memeory_for_clebsch_gordan_coefficients(const size_t bandlimit);

void free_memory_for_clebsch_gordan_coefficients(const size_t bandlimit, cg_table clebsch_gordan_table);

bool validate_clebsch_gordan_table_order(   const size_t    bandlimit, 
                                            const long int  l1, 
                                            const long int  l2, 
                                            const long int  l, 
                                            const long int  m, 
                                            const long int  m1);

double get_clebsch_gordan_coefficient(  const cg_table  clebsch_gordan_table, 
                                        const long int  l1, 
                                        const long int  l2, 
                                        const long int  l, 
                                        const long int  m, 
                                        const long int  m1);

void set_clebsch_gordan_coefficient(    cg_table        clebsch_gordan_table, 
                                        const long int  l1, 
                                        const long int  l2, 
                                        const long int  l, 
                                        const long int  m, 
                                        const long int  m1, 
                                        const double  coefficient);

long int maximum(long int a, long int b);

long int minimum(long int a, long int b);

long int absolute_value(long int a);
