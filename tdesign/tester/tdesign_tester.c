// TODO: write docs in doxygen
/* 
USE:
        tdesign_tester bandlimit
    where 
        bandlimit is non-negative integer
*/

#include <stdint.h>
#include <stdio.h>
#include "tdesign.h"

int main(int argc, char *argv[])
{
    if (argc==2)
    {
        size_t bandlimit;
        sscanf(argv[1], "%zu", &bandlimit);
        tdesign_cart td = read_tdesign(bandlimit);

        print_tdesign(&td);
        
        return 0;
    }
    else
    {
        printf("Wrong number of arguments.");
    }
}
