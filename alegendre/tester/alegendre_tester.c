// TODO: docs in doxygen

#include <stdio.h>
#include <math.h>
#include "alegendre.h"

int main(int argc, char *argv[])
{
    double out;
    printf("Begin init.\n");
    alegendre_eval_init_wrapper(&out);
    printf("End init.\n\n");

    if (argc==4)
    {
        long l, m;
        double x;
        sscanf(argv[1], "%ld", &l);
        sscanf(argv[2], "%ld", &m);
        sscanf(argv[3], "%lf", &x);

        printf("l = %d, m = %d, x = %f, acos(x) = %f, val = %f\n", l, m, x, acos(x), alegendre2(l, m, x));
    }
    else
    {
        printf("Wrong number of arguments.\n\n");
    }
}