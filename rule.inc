#Useful variable definitions
CDir=c
MEXDir=mex
COMONCOMPFLAGS=-Ofast
SIMDCOMPFLAGS=-msse2 -DHAVE_SSE2
DSFMTVAR=-DSFMT_MEXP=2281


SFMTCOMPSTR=COMPFLAGS='$COMPFLAGS $(COMONCOMPFLAGS) $(SIMDCOMPFLAGS)' $(DSFMTVAR)
NONSFMTCOMPSTR=COMPFLAGS='$COMPFLAGS $(COMONCOMPFLAGS)'

SFMTFILE=$(SMFTDir)/SFMT.c
ROTATIONFILES=$(ROTATIONSAMPLINGDir)/$(CDir)/utility_functions.c \
				$(ROTATIONSAMPLINGDir)/$(CDir)/healpix_functions.c \
				$(ROTATIONSAMPLINGDir)/$(CDir)/rotation_grid.c
CGFILE=$(CGDir)/$(CDir)/clebsch_gordan_coefficients.c
SHFILE=$(SHDir)/$(CDir)/spherical_harmonics.c
SPECTRAFILE=$(SPECTRADir)/$(CDir)/spherical_spectra.c

INCLUDEPATH=-I$(SMFTDir) -I$(ROTATIONSAMPLINGDir)/$(CDir) -I$(ALIGNMENTDir)/$(CDir) -I$(CGDir)/$(CDir) -I$(ESTIMATIONDir)/$(CDir) -I$(INVERSIONDir)/$(CDir) -I$(SHDir)/$(CDir) -I$(SPECTRADir)/$(CDir) -I$(TDESIGNDir)/$(CDir)
OUTPATH=-outdir $(MEXOUTDIR)

#Rule names
ROTATIONSAMPLINGRULES=applyRotation \
			generateUniformlyRandomRotations \
			generateUniformRotationSequencePoint
CLEBSCHGORDANRULES=cgVector
SHRULES=normalizeSHC \
			randomNormalizedSHC
SPECTRARULES=calculateBispectrum \
			calculateBispetrumGradient \
			powerSpectrum

#Rules
compile : 	$(ROTATIONSAMPLINGRULES) \
				$(CLEBSCHGORDANRULES) \
				$(SHRULES) \
				$(SPECTRARULES)


applyRotation : 
	$(MEX) $(MATLABFAGS) $(SFMTCOMPSTR) $(OUTPATH) $(INCLUDEPATH) $(ROTATIONSAMPLINGDir)/$(MEXDir)/applyRotation.c $(ROTATIONFILES) $(SFMTFILE) 

generateUniformlyRandomRotations : 
	$(MEX) $(MATLABFAGS) $(SFMTCOMPSTR) $(OUTPATH) $(INCLUDEPATH) $(ROTATIONSAMPLINGDir)/$(MEXDir)/generateUniformlyRandomRotations.c $(ROTATIONFILES) $(SFMTFILE)
	
generateUniformRotationSequencePoint : 
	$(MEX) $(MATLABFAGS) $(SFMTCOMPSTR) $(OUTPATH) $(INCLUDEPATH) $(ROTATIONSAMPLINGDir)/$(MEXDir)/generateUniformRotationSequencePoint.c $(ROTATIONFILES) $(SFMTFILE)


cgVector :
	$(MEX) $(MATLABFAGS) $(NONSFMTCOMPSTR) $(OUTPATH) $(INCLUDEPATH) $(CGDir)/$(MEXDir)/cgVector.c $(CGFILE)
	

normalizeSHC :
	$(MEX) $(MATLABFAGS) $(SFMTCOMPSTR) $(OUTPATH) $(INCLUDEPATH) $(SHDir)/$(MEXDir)/normalizeSHC.c $(SHFILE) $(SFMTFILE)

randomNormalizedSHC :
	$(MEX) $(MATLABFAGS) $(SFMTCOMPSTR) $(OUTPATH) $(INCLUDEPATH) $(SHDir)/$(MEXDir)/randomNormalizedSHC.c $(SHFILE) $(SFMTFILE) 



calculateBispectrum :
	$(MEX) $(MATLABFAGS) $(SFMTCOMPSTR) $(OUTPATH) $(INCLUDEPATH) $(SPECTRADir)/$(MEXDir)/calculateBispectrum.c $(CGFILE) $(SHFILE) $(SPECTRAFILE) $(SFMTFILE)

calculateBispetrumGradient :
	$(MEX) $(MATLABFAGS) $(SFMTCOMPSTR) $(OUTPATH) $(INCLUDEPATH) $(SPECTRADir)/$(MEXDir)/calculateBispetrumGradient.c $(CGFILE) $(SHFILE) $(SPECTRAFILE) $(SFMTFILE)

powerSpectrum :
	$(MEX) $(MATLABFAGS) $(SFMTCOMPSTR) $(OUTPATH) $(INCLUDEPATH) $(SPECTRADir)/$(MEXDir)/powerSpectrum.c $(CGFILE) $(SHFILE) $(SPECTRAFILE) $(SFMTFILE)


.PHONY : clean
clean :
	rm -r bin/*
