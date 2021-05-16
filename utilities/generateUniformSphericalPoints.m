function out = generateUniformSphericalPoints(pointNo)
% out = generateUniformSphericalPoints(pointNo)     generates pointNo
%           number of points approximately uniformly distrubted on the unit
%           sphere.
% 
% REFERENCE
%   The code is based on this answer from StackExchagne:
%       https://stackoverflow.com/a/44164075
% 

indices = 0.5:1:(pointNo-0.5);
phi = acos(1 - 2*indices/pointNo);
theta = (1 + sqrt(5)) * (pi.* indices);

x = cos(theta).*sin(phi);
x = x(:);

y = sin(theta).*sin(phi);
y = y(:);

z = cos(phi);
z = z(:);

out = [x, y, z];
