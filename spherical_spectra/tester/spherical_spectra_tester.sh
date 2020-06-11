#!/bin/bash
gcc -msse2 -DHAVE_SSE2 \
    -DSFMT_MEXP=2281 \
    -I../c \
    -I../../extern/SFMT \
    -I../../spherical_harmonics/c \
    -I../../clebsch_gordan/c \
        spherical_spectra_tester.c \
        ../../clebsch_gordan/c/clebsch_gordan_coefficients.c \
        ../../spherical_harmonics/c/spherical_harmonics.c \
        ../../spherical_spectra/c/spherical_spectra.c \
        ../../extern/SFMT/SFMT.c \
    -o spherical_spectra_tester.exe
