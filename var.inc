#Compilation variables
CC=gcc
FC=gfortran
MATLABFAGS=-R2018a
MEX=mex.bat
OUTDIRCMD=-o

#Include directory variables
SMFTDir=extern/SFMT
ROTATIONSAMPLINGDir=extern/RotationSampling
ALIGNMENTDir=alignment
CGDir=clebsch_gordan
ESTIMATIONDir=estimation
INVERSIONDir=inversion
SHDir=spherical_harmonics
SPECTRADir=spherical_spectra
TDESIGNDir=tdesign
ALEGDir=alegendre
ALEGEVALDir=extern/ALegendreEval
GFORTLIBPATH=/c/msys64/mingw64/lib/gcc/x86_64-w64-mingw32/10.1.0

#Output directory variables
MEXOUTDIR=bin/mex
TESTEROUTDIR=bin/testers