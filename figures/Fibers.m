fig = figure;

t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');


lambda = 1;

% Plane
nexttile;
hold;

scatter(0, 0, 20, 'filled', 'red');

theta = 0:0.05:2*pi;
x = cos(theta);
y = sin(theta);
for n=1:20
    plot(2*n*lambda*x, 2*n*lambda*y, '--r');
    plot((2*n+1)*lambda*x, (2*n+1)*lambda*y, ':b');
end

angles = 2*pi*(0:45:359)/360;
points = -20:20;
points = [zeros(size(points)); points];
for a=angles
    M = [cos(a), -sin(a); sin(a), cos(a)];
    rotPoints = M*points;
    plot(rotPoints(1, :), rotPoints(2, :), 'm');
end

hold off;

xlim([-10, 10]);
ylim([-10, 10]);

xticks([]);
yticks([]);

set(gca, 'XAxisLocation', 'origin');
set(gca, 'YAxisLocation', 'origin');

title('Plane');

% Sphere
nexttile;
hold;

[x, y, z] = sphere(128);
s = surfl(x, y, z); 
colormap('gray');
set(s, 'FaceAlpha', 0.25)
shading FLAT;

np = [0, 0, 1];
scatter3(np(1), np(2), np(3), 30, 'filled', 'r');
scatter3(-np(1), -np(2), -np(3), 30, 'filled', 'b');

theta = 0:0.05:2*pi;
greatCircle = [cos(theta); sin(theta)];
greatCircle = [zeros(1, size(greatCircle, 2)); greatCircle];
angles = 2*pi*(0:45:359)/360;
for a=angles
    M = eye(3);
    M(1:2, 1:2) = [cos(a), -sin(a); sin(a), cos(a)];
    rotGreatCircle = M*greatCircle;
    plot3(rotGreatCircle(1, :), rotGreatCircle(2, :), rotGreatCircle(3, :), 'm');
end


hold off;

axis off;
view(5, 50)
xlim([-1, 1]);
ylim([-1, 1]);

title('Sphere');
