#Compilation variables
CC=gcc
FC=gfortran
MATLABFLAGS=-R2018a CFLAGS='$$CFLAGS -Ofast'
MATLABFLAGSBISP=-R2018a CFLAGS='$$CFLAGS -Ofast --param early-inlining-insns=25 --param max-inline-insns-auto=15'
#-fopt-info-missed 
MEX=mex.bat
