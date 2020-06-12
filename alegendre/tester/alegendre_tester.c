// TODO

#include <stdio.h>

void alegendre_eval_init_wrapper(double *dsize);

int main()
{
    double out;
    printf("Begin init.\n");
    alegendre_eval_init_wrapper(&out);
    printf("End init.\n");
}