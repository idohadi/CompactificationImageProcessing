include make.inc

EXTERNAL = 	FastSphericalHarmonicsTransform \
			SmallRotationToolbox

all : EXTERNAL

# External repos
FastSphericalHarmonicsTransform :
	$(make) -C extern/FastSphericalHarmonicsTransform
	
SmallRotationToolbox :
	$(make) -C extern/SmallRotationToolbox
	
