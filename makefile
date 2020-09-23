include make.inc

EXTERNAL = 	SmallRotationToolbox

MEXFILES = 	normSHC_mex \
			ClebschGordanCoeffs_mex \
			buildCGTable_mex \
			bispectrum_mex \
			powerSpectrum_mex
#TODO: add the estimation mex funcs once I update their Clebsch-Gordan code base
			
all : EXTERNAL MEXFILES

# External repos	
SmallRotationToolbox :
	$(make) -C extern/SmallRotationToolbox
	

# MEX files
normSHC_mex : 
	$(MEX) $(MATLABFLAGS) spherical_harmonics/normSHC_mex.c -outdir spherical_harmonics/

ClebschGordanCoeffs_mex : 
	$(MEX) $(MATLABFLAGS) \
		clebsch_gordan/ClebschGordanCoeffs_mex.c \
		clebsch_gordan/clebsch_gordan_coefficients.c \
		-outdir clebsch_gordan/

buildCGTable_mex : 
	$(MEX) $(MATLABFLAGS) \
		clebsch_gordan/buildCGTable_mex.c \
		clebsch_gordan/clebsch_gordan_coefficients.c \
		-outdir clebsch_gordan/

powerSpectrum_mex : 
	$(MEX) $(MATLABFLAGS) spherical_spectra/powerSpectrum_mex.c -outdir spherical_spectra/

bispectrum_mex : 
	$(MEX) $(MATLABFLAGSBISP) spherical_spectra/bispectrum_mex.c -outdir spherical_spectra/

buildV_mex : 
	$(MEX) $(MATLABFLAGS) estimation/buildV_mex.c -outdir estimation/
	
buildK_mex : 
	$(MEX) $(MATLABFLAGS) estimation/buildK_mex.c -outdir estimation/
	

# Cleaning rules
.PHONY : clean
clean : clean_o clean_mex clean_external
	
clean_o :
	rm -f *.o
	
clean_mex :
	rm -f *.mexw64
	
clean_external :
	$(MAKE) clean -C $(CURDIR)/extern/SmallRotationToolbox
	