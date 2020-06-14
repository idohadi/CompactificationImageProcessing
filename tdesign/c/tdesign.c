// TODO: write docs in doxygen

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "tdesign.h"

// TODO: add list of file names

// A list of t-design file names
char * const bandlimit_filename[]
    = 
    {
        "../data/sf001.00003",
        "../data/sf002.00006",
        "../data/sf003.00008",
        "../data/sf004.00014",
        "../data/sf005.00018",
        "../data/sf006.00026",
        "../data/sf007.00032",
        "../data/sf008.00042",
        "../data/sf009.00050",
        "../data/sf010.00062",
        "../data/sf011.00072",
        "../data/sf012.00086",
        "../data/sf013.00098",
        "../data/sf014.00114",
        "../data/sf015.00128",
        "../data/sf016.00146",
        "../data/sf017.00163",
        "../data/sf018.00182",
        "../data/sf019.00201",
        "../data/sf020.00222",
        "../data/sf021.00243",
        "../data/sf022.00266",
        "../data/sf023.00289",
        "../data/sf024.00314",
        "../data/sf025.00339",
        "../data/sf026.00366",
        "../data/sf027.00393",
        "../data/sf028.00422",
        "../data/sf029.00451",
        "../data/sf030.00482",
        "../data/sf031.00513",
        "../data/sf032.00546",
        "../data/sf033.00579",
        "../data/sf034.00614",
        "../data/sf035.00649",
        "../data/sf036.00686",
        "../data/sf037.00723",
        "../data/sf038.00762",
        "../data/sf039.00801",
        "../data/sf040.00842",
        "../data/sf041.00883",
        "../data/sf042.00926",
        "../data/sf043.00969",
        "../data/sf044.01014",
        "../data/sf045.01059",
        "../data/sf046.01106",
        "../data/sf047.01153",
        "../data/sf048.01202",
        "../data/sf049.01251",
        "../data/sf050.01302",
        "../data/sf051.01353",
        "../data/sf052.01406",
        "../data/sf053.01459",
        "../data/sf054.01514",
        "../data/sf055.01569",
        "../data/sf056.01626",
        "../data/sf057.01683",
        "../data/sf058.01742",
        "../data/sf059.01801",
        "../data/sf060.01862",
        "../data/sf061.01923",
        "../data/sf062.01986",
        "../data/sf063.02049",
        "../data/sf064.02114",
        "../data/sf065.02179",
        "../data/sf066.02246",
        "../data/sf067.02313",
        "../data/sf068.02382",
        "../data/sf069.02451",
        "../data/sf070.02522",
        "../data/sf071.02593",
        "../data/sf072.02666",
        "../data/sf073.02739",
        "../data/sf074.02814",
        "../data/sf075.02889",
        "../data/sf076.02966",
        "../data/sf077.03043",
        "../data/sf078.03122",
        "../data/sf079.03201",
        "../data/sf080.03282",
        "../data/sf081.03363",
        "../data/sf082.03446",
        "../data/sf083.03529",
        "../data/sf084.03614",
        "../data/sf085.03699",
        "../data/sf086.03786",
        "../data/sf087.03873",
        "../data/sf088.03962",
        "../data/sf089.04051",
        "../data/sf090.04142",
        "../data/sf091.04233",
        "../data/sf092.04326",
        "../data/sf093.04419",
        "../data/sf094.04514",
        "../data/sf095.04609",
        "../data/sf096.04706",
        "../data/sf097.04803",
        "../data/sf098.04902",
        "../data/sf099.05001",
        "../data/sf100.05102",
        "../data/sf101.05203",
        "../data/sf102.05306",
        "../data/sf103.05409",
        "../data/sf104.05514",
        "../data/sf105.05619",
        "../data/sf106.05726",
        "../data/sf107.05833",
        "../data/sf108.05942",
        "../data/sf109.06051",
        "../data/sf110.06162",
        "../data/sf111.06273",
        "../data/sf112.06386",
        "../data/sf113.06499",
        "../data/sf114.06614",
        "../data/sf115.06729",
        "../data/sf116.06846",
        "../data/sf117.06963",
        "../data/sf118.07082",
        "../data/sf119.07201",
        "../data/sf120.07322",
        "../data/sf121.07443",
        "../data/sf122.07566",
        "../data/sf123.07689",
        "../data/sf124.07814",
        "../data/sf125.07939",
        "../data/sf126.08066",
        "../data/sf127.08193",
        "../data/sf128.08322",
        "../data/sf129.08451",
        "../data/sf130.08582",
        "../data/sf131.08713",
        "../data/sf132.08846",
        "../data/sf133.08979",
        "../data/sf134.09114",
        "../data/sf135.09249",
        "../data/sf136.09386",
        "../data/sf137.09523",
        "../data/sf138.09662",
        "../data/sf139.09801",
        "../data/sf140.09942",
        "../data/sf141.10083",
        "../data/sf142.10226",
        "../data/sf143.10369",
        "../data/sf144.10514",
        "../data/sf145.10659",
        "../data/sf146.10806",
        "../data/sf147.10953",
        "../data/sf148.11102",
        "../data/sf149.11251",
        "../data/sf150.11402",
        "../data/sf151.11553",
        "../data/sf152.11706",
        "../data/sf153.11859",
        "../data/sf154.12014",
        "../data/sf155.12169",
        "../data/sf156.12326",
        "../data/sf157.12483",
        "../data/sf158.12642",
        "../data/sf159.12801",
        "../data/sf160.12962",
        "../data/sf161.13123",
        "../data/sf162.13286",
        "../data/sf163.13449",
        "../data/sf164.13614",
        "../data/sf165.13779",
        "../data/sf166.13946",
        "../data/sf167.14113",
        "../data/sf168.14282",
        "../data/sf169.14451",
        "../data/sf170.14622",
        "../data/sf171.14793",
        "../data/sf172.14966",
        "../data/sf173.15139",
        "../data/sf174.15314",
        "../data/sf175.15489",
        "../data/sf176.15666",
        "../data/sf177.15843",
        "../data/sf178.16022",
        "../data/sf179.16201",
        "../data/sf180.16382"
    };

long const tdesign_length[]
    = 
    {
        3,
        6,
        8,
        14,
        18,
        26,
        32,
        42,
        50,
        62,
        72,
        86,
        98,
        114,
        128,
        146,
        163,
        182,
        201,
        222,
        243,
        266,
        289,
        314,
        339,
        366,
        393,
        422,
        451,
        482,
        513,
        546,
        579,
        614,
        649,
        686,
        723,
        762,
        801,
        842,
        883,
        926,
        969,
        1014,
        1059,
        1106,
        1153,
        1202,
        1251,
        1302,
        1353,
        1406,
        1459,
        1514,
        1569,
        1626,
        1683,
        1742,
        1801,
        1862,
        1923,
        1986,
        2049,
        2114,
        2179,
        2246,
        2313,
        2382,
        2451,
        2522,
        2593,
        2666,
        2739,
        2814,
        2889,
        2966,
        3043,
        3122,
        3201,
        3282,
        3363,
        3446,
        3529,
        3614,
        3699,
        3786,
        3873,
        3962,
        4051,
        4142,
        4233,
        4326,
        4419,
        4514,
        4609,
        4706,
        4803,
        4902,
        5001,
        5102,
        5203,
        5306,
        5409,
        5514,
        5619,
        5726,
        5833,
        5942,
        6051,
        6162,
        6273,
        6386,
        6499,
        6614,
        6729,
        6846,
        6963,
        7082,
        7201,
        7322,
        7443,
        7566,
        7689,
        7814,
        7939,
        8066,
        8193,
        8322,
        8451,
        8582,
        8713,
        8846,
        8979,
        9114,
        9249,
        9386,
        9523,
        9662,
        9801,
        9942,
        10083,
        10226,
        10369,
        10514,
        10659,
        10806,
        10953,
        11102,
        11251,
        11402,
        11553,
        11706,
        11859,
        12014,
        12169,
        12326,
        12483,
        12642,
        12801,
        12962,
        13123,
        13286,
        13449,
        13614,
        13779,
        13946,
        14113,
        14282,
        14451,
        14622,
        14793,
        14966,
        15139,
        15314,
        15489,
        15666,
        15843,
        16022,
        16201,
        16382
    };

char *bandlimit_to_filename(const size_t bandlimit)
{
    if (bandlimit>=1 && bandlimit<=180)
    {
        return bandlimit_filename[bandlimit-1];
    }
    return NULL;
}

long bandlimit_to_tdesign_length(const size_t bandlimit)
{
    if (bandlimit>=1 && bandlimit<=180)
    {
        return tdesign_length[bandlimit-1];
    }
    return 0;
}

double *allocate_tdesign(const size_t bandlimit, COORD_SYSTEM sys)
{
    return malloc((sys==CART ? 3 : 2)*bandlimit_to_tdesign_length(bandlimit)*sizeof(double));
}

void deallocate_tdesign(double *tdesign)
{
    free(tdesign);
}

tdesign_cart read_tdesign(const size_t bandlimit)
{
    FILE *f = fopen(bandlimit_to_filename(bandlimit), "r");
    if (f!=NULL)
    {
        tdesign_cart td;
        td.bandlimit = bandlimit;
        td.length = bandlimit_to_tdesign_length(bandlimit);
        td.tdesign = allocate_tdesign(bandlimit, CART);

        double x1, x2, x3;
        size_t row = 0;
        while (fscanf(f, " %23lf %23lf %23lf", &x1, &x2, &x3)!=EOF)
        {
            td.tdesign[3*row] = x1;
            td.tdesign[3*row+1] = x2;
            td.tdesign[3*row+2] = x3;
            ++row;
        }
        return td;
    }
    else
    {
        printf("Reading t-design unsuccessful.\n");
    }
}


void print_tdesign(tdesign_cart *td)
{
    printf("Bandlimit = %d. t-design length = %d.\n", td->bandlimit, td->length);
    printf("=========================================================================\n");
    printf("\tRow\t\t   x1\t\t   x2\t\t   x3\n");
    printf("=========================================================================\n");
    for (long row = 0; row<td->length; ++row)
    {
        printf("\t%d\t\t% 10.10lf\t% 10.10lf\t% 10.10lf\n", row+1, td->tdesign[3*row], td->tdesign[3*row+1], td->tdesign[3*row+2]);
    }
}
