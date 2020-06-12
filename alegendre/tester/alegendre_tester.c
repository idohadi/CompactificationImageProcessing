// TODO

#include <stdio.h>
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
        double t;
        sscanf(argv[1], "%ld", &l);
        sscanf(argv[2], "%ld", &m);
        sscanf(argv[3], "%lf", &t);

        printf("l = %d, m = %d, t = %f, val = %f\n", l, m, t, alegendre(l, m, t));
    }
    else
    {
        printf("Wrong number of arguments.\n\n");
    }
}