
const char *bandlimit_filename[];
const long tdesign_length[];

/* t-design in Cartesian coordinates. */
typedef 
    struct tdesign_cart
    {
        size_t bandlimit;
        double *tdesign;
    }
    tdesign_cart;

/* t-design in spherical coordinates. */
typedef 
    struct tdesign_sph
    {
        size_t bandlimit;
        double *tdesign;
    }
    tdesign_sph;

typedef enum COORD_SYSTEM {EUC, SPH} COORD_SYSTEM;


/* Functions */
char *bandlimit_to_filename(const size_t bandlimit);

long bandlimit_to_tdesign_length(const size_t bandlimit);

double *allocate_tdesign(const size_t bandlimit, COORD_SYSTEM sys);

void deallocate_tdesign(double *tdesign);
