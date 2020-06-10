// TODO

#include <stdint.h>

/*  
    for real-valued spherical functions: spherical harmonics coefficients format is
        f_{r, 0, 0}, 
            f_{r, 1, 0}, f_{r, 1, 1}, f_{i, 1, 1}, 
                f_{r, 2, 0}, f_{r, 2, 1}, f_{i, 2, 1}, f_{r, 2, 2}, f_{i, 2, 2}, ...
                    f_{r, bandlimit, 0}, f_{r, bandlimit, 1}, f_{i, bandlimit, 1}, ..., f_{r, bandlimit, bandlimit}, f_{i, bandlimit, bandlimit}
*/
typedef 
    struct r_shc
    {
        double *coefficients;
        size_t bandlimit;
    }
    r_shc;

/*  
    for complex-valued spherical functions: spherical harmonics coefficients format is
        f_{r, 0, 0}, 
            f_{r, 1, -1}, f_{i, 1, -1}, f_{r, 1, 0}, f_{i, 1, 0}, f_{r, 1, 1}, f_{i, 1, 1}, 
                f_{r, 2, -2}, f_{i, 2, -2}, f_{r, 2, -1}, f_{i, 2, -1}, f_{r, 2, 0}, f_{i, 2, 0}, f_{r, 2, 1}, f_{i, 2, 1}, f_{r, 2, 2}, f_{i, 2, 2}, 
                ...
*/
typedef 
    struct c_shc
    {
        double *coefficients;
        size_t bandlimit;
    }
    c_shc;

typedef enum PART {REAL_PART, IMAG_PART} PART;
