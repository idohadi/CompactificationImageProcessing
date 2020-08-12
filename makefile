include make.inc

EXTERNAL = 	FastSphericalHarmonicsTransform \
			SmallRotationToolbox

MEXFILES = 	normSHC_mex \
			ClebschGordanCoeffs_mex

all : EXTERNAL MEXFILES

# External repos
FastSphericalHarmonicsTransform :
	$(make) -C extern/FastSphericalHarmonicsTransform
	
SmallRotationToolbox :
	$(make) -C extern/SmallRotationToolbox
	

# MEX files
normSHC_mex : 
	$(MEX) $(MATLABFAGS) spherical_harmonics/normSHC_mex.c -outdir spherical_harmonics/

ClebschGordanCoeffs_mex : 
	$(MEX) $(MATLABFAGS) \
		clebsch_gordan/ClebschGordanCoeffs_mex.c \
		clebsch_gordan/clebsch_gordan_coefficients.c \
		-outdir clebsch_gordan/

buildCGTable_mex : 
	$(MEX) $(MATLABFAGS) \
		clebsch_gordan/buildCGTable_mex.c \
		clebsch_gordan/clebsch_gordan_coefficients.c \
		-outdir clebsch_gordan/

powerSpectrum_mex : 
	$(MEX) $(MATLABFAGS) spherical_spectra/powerSpectrum_mex.c -outdir spherical_spectra/

bispectrum_mex : 
	$(MEX) $(MATLABFAGS) spherical_spectra/bispectrum_mex.c -outdir spherical_spectra/
	

# Cleaning rules
.PHONY : clean
clean : clean_mex clean_external
	
clean_mex :
	rm -f *.mexw64
	
clean_external :
	$(MAKE) clean -C $(CURDIR)/extern/SmallRotationToolbox
	$(MAKE) clean -C $(CURDIR)/extern/FastSphericalHarmonicsTransform
