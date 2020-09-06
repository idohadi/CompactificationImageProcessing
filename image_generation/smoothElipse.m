function im = smoothElipse(a, b, imageSize, bandlimit)

%% Input validation
narginchk(2, 4);

if nargin==2
    imageSize = 151;
    bandlimit = 16;
elseif nargin==3
    bandlimit = 16;
end

assert(round(imageSize)==imageSize ...
    && imageSize>0 ...
    && numel(imageSize)==1, ...
    'Image size must be a positive integer.');
assert(round(bandlimit)==bandlimit ...
    && bandlimit>0 ...
    && numel(bandlimit)==1, ...
    'Bandlimit size must be a positive integer scalar.');

%% Generate ellipse
im = zeros(imageSize, imageSize);

Sp = linspace(-0.5, 0.5, imageSize);
[X, Y] = meshgrid(Sp, Sp);

mask = X.^2/a + Y.^2/b <= 0.5;

im(mask) = 1;
im = imgaussfilt(im, 5);

td = loadtd('sf045.01059');
shc = image2shc(im, 16, td, KondorBackProj(1), [-0.5, 0.5], [-0.5, 0.5]);
im = shc2image(shc, 16, imageSize, KondorProj(1), [-0.5, 0.5], [-0.5, 0.5]);
im = real(im);
