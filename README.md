# Compactification of the Rigid Motions Group in Image Processing

Code repository of the paper [Compactification of the Rigid Motions Group in Image Processing](https://arxiv.org/abs/2106.13505).

# Installation

**Note.** The code depends on a submodule [SmallRotationToolbox](https://github.com/idohadi/SmallRotationToolbox), which is assumed to be placed in the __extern__ folder. 

1. In order to clone both this repository and its submodule, run 

    git clone --recurse-submodules https://github.com/idohadi/CompactificationImageProcessing.git
    
1. Change folder to root of the repository and run `make`, to compile all MATLAB's mex function.

1. Run __setup.m__ to add all necessary folders to MATLAB's path.

**Note.** MEX function compilation was tested only with __gcc__.

# Clebsch-Gordan Coefficients
In order to precompute the Clebsch-Gordan coefficients, create the folder __clebsch_gordan/ClebschGordanCoeffs/__ and run `buildCGTable(L)`, where `L` is the bandlimit of the built table.

**Note.** In our numerical experiments we used `L=16` or `L=30`.
