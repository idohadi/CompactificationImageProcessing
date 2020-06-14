#Useful variable definitions
CDir=c
FDir=f
MEXDir=mex
TESTDir=tester
COMONCOMPFLAGS=-Ofast
FCOMONCOMPFLAGS=-Ofast
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
TDFILE=$(TDESIGNDir)/$(CDir)/tdesign.c

INCLUDEPATH=-I$(SMFTDir) \
	-I$(ROTATIONSAMPLINGDir)/$(CDir) \
	-I$(ALIGNMENTDir)/$(CDir) \
	-I$(CGDir)/$(CDir) \
	-I$(ESTIMATIONDir)/$(CDir) \
	-I$(SHDir)/$(CDir) \
	-I$(SPECTRADir)/$(CDir) \
	-I$(TDESIGNDir)/$(CDir) \
	-I$(ALEGDir)/$(FDir) -I$(ALEGDir)/$(CDir)
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
ALEGENDRERULES=alegendre_mex
TESTERSRULES=clebsch_gordan_tester \
			spherical_harmonics_tester \
			spherical_spectra_tester \
			alegendre_tester

#Rules
compile : 	$(ROTATIONSAMPLINGRULES) \
				$(CLEBSCHGORDANRULES) \
				$(SHRULES) \
				$(SPECTRARULES) \
				$(ALEGENDRERULES) \
				$(TESTERSRULES)


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



alegendre_mex :
	mkdir -p $(MEXOUTDIR)
	$(FC) -c $(INCLUDEPATH) $(FCOMONCOMPFLAGS) $(ALEGEVALDir)/bessel_eval.f90 -o $(MEXOUTDIR)/bessel_eval.o -J$(MEXOUTDIR)
	$(FC) -c $(INCLUDEPATH) $(FCOMONCOMPFLAGS) $(ALEGEVALDir)/alegendre_eval.f90 -o $(MEXOUTDIR)/alegendre_eval.o -J$(MEXOUTDIR)
	$(FC) -c $(INCLUDEPATH) $(FCOMONCOMPFLAGS) $(ALEGDir)/$(FDir)/alegendre_wrapper.f90 -o $(MEXOUTDIR)/alegendre_wrapper.o -J$(MEXOUTDIR)
	$(CC) -c $(INCLUDEPATH) $(COMONCOMPFLAGS) $(ALEGDir)/$(CDir)/alegendre.c -o $(MEXOUTDIR)/alegendre.o -J$(MEXOUTDIR)
	
	$(MEX) $(MATLABFAGS) $(NONSFMTCOMPSTR) $(OUTPATH) $(INCLUDEPATH) $(ALEGDir)/$(MEXDir)/evalAlegendre.c \
		$(MEXOUTDIR)/alegendre_wrapper.o \
		$(MEXOUTDIR)/alegendre.o \
		$(MEXOUTDIR)/alegendre_eval.o \
		$(MEXOUTDIR)/bessel_eval.o \
		-L$(GFORTLIBPATH) -lgfortran -lquadmath
	
	rm $(MEXOUTDIR)/*.mod
	rm $(MEXOUTDIR)/*.o
	cp $(ALEGEVALDir)/*.bin $(MEXOUTDIR)
	cp $(ALEGEVALDir)/*.bin1 $(MEXOUTDIR)
	cp $(ALEGEVALDir)/*.bin2 $(MEXOUTDIR)
	cp $(ALEGEVALDir)/*.bin3 $(MEXOUTDIR)


.PHONY3 : create_testers_folder
create_testers_folder : 
	mkdir -p $(TESTEROUTDIR)

.PHONY2 : testers
testers : $(TESTERSRULES)
	
clebsch_gordan_tester : create_testers_folder
	$(CC) $(INCLUDEPATH) $(CGDir)/$(TESTDir)/clebsch_gordan_tester.c $(CGFILE) $(OUTDIRCMD) $(TESTEROUTDIR)/clebsch_gordan_tester.exe

spherical_harmonics_tester : create_testers_folder
	$(CC) $(INCLUDEPATH) $(SHDir)/$(TESTDir)/spherical_harmonics_tester.c $(CGFILE) $(SHFILE) $(SFMTFILE)  $(OUTDIRCMD) $(TESTEROUTDIR)/spherical_harmonics_tester.exe $(DSFMTVAR) $(SIMDCOMPFLAGS)

spherical_spectra_tester : create_testers_folder
	$(CC) $(INCLUDEPATH) $(SPECTRADir)/$(TESTDir)/spherical_spectra_tester.c $(CGFILE) $(SHFILE) $(SPECTRAFILE) $(SFMTFILE) $(OUTDIRCMD) $(TESTEROUTDIR)/spherical_spectra_tester.exe $(DSFMTVAR) $(SIMDCOMPFLAGS)

alegendre_tester : create_testers_folder
	$(FC) -c $(INCLUDEPATH) $(FCOMONCOMPFLAGS) $(ALEGEVALDir)/bessel_eval.f90 -o $(TESTEROUTDIR)/bessel_eval.o -J$(TESTEROUTDIR)
	$(FC) -c $(INCLUDEPATH) $(FCOMONCOMPFLAGS) $(ALEGEVALDir)/alegendre_eval.f90 -o $(TESTEROUTDIR)/alegendre_eval.o -J$(TESTEROUTDIR)
	$(FC) -c $(INCLUDEPATH) $(FCOMONCOMPFLAGS) $(ALEGDir)/$(FDir)/alegendre_wrapper.f90 -o $(TESTEROUTDIR)/alegendre_wrapper.o -J$(TESTEROUTDIR)
	$(CC) -c $(INCLUDEPATH) $(COMONCOMPFLAGS) $(ALEGDir)/$(CDir)/alegendre.c -o $(TESTEROUTDIR)/alegendre.o -J$(TESTEROUTDIR)
	$(CC) $(INCLUDEPATH) $(COMONCOMPFLAGS) $(ALEGDir)/$(TESTDir)/alegendre_tester.c \
		$(TESTEROUTDIR)/alegendre_wrapper.o \
		$(TESTEROUTDIR)/alegendre.o \
		$(TESTEROUTDIR)/alegendre_eval.o \
		$(TESTEROUTDIR)/bessel_eval.o \
		-lgfortran \
		$(OUTDIRCMD) $(TESTEROUTDIR)/alegendre_tester.exe
	rm $(TESTEROUTDIR)/*.mod
	rm $(TESTEROUTDIR)/*.o
	cp $(ALEGEVALDir)/*.bin $(TESTEROUTDIR)
	cp $(ALEGEVALDir)/*.bin1 $(TESTEROUTDIR)
	cp $(ALEGEVALDir)/*.bin2 $(TESTEROUTDIR)
	cp $(ALEGEVALDir)/*.bin3 $(TESTEROUTDIR)

tdesign_tester : create_testers_folder
	$(CC) $(INCLUDEPATH) $(TDESIGNDir)/$(TESTDir)/tdesign_tester.c $(TDFILE) $(OUTDIRCMD) $(TESTEROUTDIR)/tdesign_tester.exe
	mkdir -p bin/data
	cp $(TDESIGNDir)/data/* bin/data
	
.PHONY : clean
clean :
	rm -r bin


