// TODO

#include <stdint.h>

typedef 
    struct r_shc
    {
        double *coefficients;
        size_t bandlimit;
    }
    r_shc;

typedef 
    struct c_shc
    {
        double *coefficients;
        size_t bandlimit;
    }
    c_shc;

typedef enum PART {REAL_PART, IMAG_PART} PART;
