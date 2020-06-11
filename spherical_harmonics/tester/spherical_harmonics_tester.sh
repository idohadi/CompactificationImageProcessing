#!/bin/bash
#SFMTDir="../../extern/SFMT"
#SFMTFILE="$(SMFTDir)/SFMT.c"
gcc -msse2 -DHAVE_SSE2 \
    -DSFMT_MEXP=2281 \
    -I../c -I../../extern/SFMT \
    spherical_harmonics_tester.c ../c/spherical_harmonics.c ../../extern/SFMT/SFMT.c \
    -o spherical_harmonics_tester.exe
