function P = interp2linop(imageSize, tDesign, interval, scalingParam)
%%
% Call format
%   P = interp2linop(imageSize, tDesign, interval, scalingParam)
% 
% Compute the linear operator P satisfying 
%   P I = 2dInterpolation(KondorBackProjection(tDesign) points in the cube inerval^2)
% where I is an image.
% 
% Input arguments
%   imageSize       double      positive integer, the image is of size 
%                                   imageSize x imageSize.
%   tDesign         double      N x 3 array, a spherical design (a 
%                               t-design) in Cartesian coordiantes.
%   interval        double      1 x 2 array, interval(1) is the lower bound
%                                of the interval and interval(2).
%   scalingParam    double      positive number, scaling parameter for 
%                               the projection.
% 
% Output arguments
%   P               double      a matrix represneting the linear opeartor
%                               described above.
% 
% Notes
%   This function performs no input checks.
% 
% Reference
%   None
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2020
% ***********************************************************

%% Compute the linear opeartor
% Back-project t-design points onto R^2
[tDesignTheta, tDesignPhi, ~] ...
    = cart2sph2(tDesign(:, 1), tDesign(:, 2), tDesign(:, 3));
[R2phi, R2rho] = KondorBackProj(tDesignTheta, tDesignPhi, scalingParam);
[R2x, R2y] = pol2cart2(R2phi, R2rho);

% Isolating t-design points falling within interval^2
tDesignInCube = R2x>=interval(1) & R2x<=interval(2) ...
    & R2y>=interval(1) & R2y<=interval(2);
R2x = R2x(tDesignInCube);
R2y = R2y(tDesignInCube);

% Creating a grid
X = linspace(interval{:}, imageSize); 
Y = linspace(interval{:}, imageSize);
[X, Y] = meshgrid(X, Y);

% Value for extrapulation outside interval^2
extval = 0;

% Compute the interpolation operator
P = zeros(length(R2x), imageSize^2);
for w=1:imageSize
    S = zeros(length(R2x), imageSize);
    parfor h=1:imageSize
        imageDelta = zeros(imageSize, imageSize);
        imageDelta(h, w) = 1;
        S(:, h) = interp2(X, Y, imageDelta, R2x, R2y, 'cubic', extval);
    end
    P(:, (w-1)*imageSize+1:w*imageSize) = S;
    disp(['Completed column ', num2str(w), ' of ', num2str(imageSize), '.']);
end
