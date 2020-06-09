// TODO

#include <stdint.h>
#include "tdesign.h"

// TODO: add list of file names
const char **bandlimit_filename = {};

char *bandlimit_to_filename(size_t bandlimit)
{
    if (bandlimit>=1)
    {
        return bandlimit_filename[bandlimit-1];
    }
    return NULL;
}

void read_tdesign()
{
    // TODO
    // Read t-design from file, based on bandlimit, not file name
}
