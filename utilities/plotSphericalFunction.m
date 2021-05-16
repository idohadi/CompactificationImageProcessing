function plotSphericalFunction(funHandle, pointNo, varargin)
% PLOTSPHERICALFUNCTION(funHandle, pointNo) creates a heatmap-like plot 
%   of funHandle on the sphere.
%   
%   Call format of funHandle must be funHandle(theta, phi), returning 
%       out(k) = function(theta(k), phi(k)).
%   out(k) must be real.
% 
%   I use the convention 0<=theta<=2pi and 0<=phi<=pi.
%   
%   ponitNo is the number of uniformly distributed points on the sphere
%   that the plot will be based on.
%   
% 
% Optional input arguments (Name, Value pairs):
%   Colormap
%       void: any MATLAB colormap to use for patch face colors
%       Default: 'hot'
%   
%   Colorbar
%       bool: show a colorbar on figure?
%       Default: true
% 
%   EdgeColor
%       void: value for EdgeColor property of patch
%       Default: 'none'
% 
%   CAxis
%       1 x 2 array: input for caxis function
%       Default: nan


%% Optional input arguments handle
p = inputParser();
addParameter(p, 'Colormap', 'hot');
addParameter(p, 'Colorbar', true);
addParameter(p, 'EdgeColor', 'none');
addParameter(p, 'CAxis', nan);
parse(p, varargin{:});

chosenColormap = p.Results.Colormap;
showColorbar = p.Results.Colorbar;
chosenEdgeColor = p.Results.EdgeColor;
CAxis = p.Results.CAxis;

%% Build the plot parameters
points = generateUniformSphericalPoints(pointNo);
triangulationIndices = convhull(points);

% Convert points to spherical coordinates


% Calculate the midpoint of every triangle in the triangulation
midpoints = 1/3*(points(triangulationIndices(:, 1), :) ...
                + points(triangulationIndices(:, 2), :) ...
                + points(triangulationIndices(:, 3), :));
midpoints = midpoints./vecnorm(midpoints, 2, 2);
[theta_middle, phi_middle] = cartesianToSphericalCoordinates(midpoints);

% Estimate function on the given points
funValues = funHandle(theta_middle, phi_middle);

%% Plot the data
% Choose the colormap
colormap(chosenColormap);

% Create the patches object
patch('Vertices', points, ...
    'Faces', triangulationIndices, ...
    'FaceVertexCData', funValues, ...
    'CDataMapping', 'scaled', ...
    'FaceColor', 'flat', ...
    'EdgeColor', chosenEdgeColor);
if showColorbar
    colorbar;
end
if ~isnan(CAxis)
    caxis(CAxis);
end
view(3);


%% Converts Cartesian coordinates into spherical coordinates
function [theta, phi] = cartesianToSphericalCoordinates(cartesianPoints)
% Given N x 3 array cartesianPoints, calcualte their representation in
% spherical coordinates.
% 
% I use the convention 0<=theta<=2pi and 0<=phi<=pi.
% 

% Convert points to spherical coordinates
x = cartesianPoints(:, 1);
y = cartesianPoints(:, 2);
z = cartesianPoints(:, 3);

theta = acos(z);

x1_proj = zeros(length(x), 1);
x2_proj = zeros(length(y), 1);

temp = sin(theta);
x1_proj(temp~=0) = x(temp~=0)./temp(temp~=0);
x2_proj(temp~=0) = y(temp~=0)./temp(temp~=0);
clear temp;

phi = atan2(x2_proj, x1_proj);
% Fix MATLAB's atan2 image being [-pi, pi]
phi(phi<0) = 2*pi + phi(phi<0);

