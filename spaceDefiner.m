spaceDefiner

% create the surface of a cylinder
Y = 0:1:200;
% Y=-0.32*(t).^2+200;
X = ((Y-200)/-0.32).^(1/2);
figure(31)
scatter(X,Y)
[X,Y,Z] = cylinder(X,20);
figure(30)
surf(X,Y,Z)
figure(31)
mesh(X,Y,Z)

% reshape the matrix into n*1 dimentions
X = reshape(X,[numel(X),1])
Y = reshape(Y,[numel(Y),1])
Z = reshape(Z,[numel(Z),1])
Z = Z*200
figure(32)
scatter3(X,Y,Z)

% Generate random points
n = 10000; % the total number of random points
xy = rand(2, n).*50-25; % Generate random points at xyc

z = rand(1, n).*200; % Generate random points at z
tri = delaunayn([X,Y,Z]); % Generate delaunay triangulization
tn = tsearchn([X,Y,Z], tri, [xy;z]'); % Determine which triangle point is within
IsInside = ~isnan(tn) % Convert to logical vector

% Only retain the points within the cylinder
xy = xy(:, IsInside);
z = z(:, IsInside);
figure(33)
scatter3(xy(1,:), xy(2,:), z)
scatter(xy(1,:), xy(2,:))