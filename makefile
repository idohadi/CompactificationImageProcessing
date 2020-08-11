include make.inc

EXTERNAL = 	FastSphericalHarmonicsTransform \
			SmallRotationToolbox

MEXFILES = 	normSHC_mex 

all : EXTERNAL

# External repos
FastSphericalHarmonicsTransform :
	$(make) -C extern/FastSphericalHarmonicsTransform
	
SmallRotationToolbox :
	$(make) -C extern/SmallRotationToolbox
	

# MEX files
normSHC_mex : 
	$(MEX) $(MATLABFAGS) spherical_harmonics/normSHC_mex.c -outdir spherical_harmonics/


# Cleaning rules
.PHONY : clean
clean : clean_mex clean_external
	
clean_mex :
	rm -f *.mexw64
	
clean_external :
	$(MAKE) clean -C $(CURDIR)/extern/SmallRotationToolbox
	$(MAKE) clean -C $(CURDIR)/extern/FastSphericalHarmonicsTransform
