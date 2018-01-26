% create the surface of a cylinder
t = 0:pi/10:pi;
[X,Y,Z] = cylinder(pi/4+cos(t));
figure(30)
surf(X,Y,Z)
figure(31)
mesh(X,Y,Z)

% reshape the matrix into n*1 dimentions
X = reshape(X,[numel(X),1])
Y = reshape(Y,[numel(Y),1])
Z = reshape(Z,[numel(Z),1])
figure(32)
scatter3(X,Y,Z)

% Generate random points
n = 1000000; % the total number of random points
xy = rand(2, n).*4-2; % Generate random points at xy
z = rand(1, n); % Generate random points at z
tri = delaunayn([X,Y,Z]); % Generate delaunay triangulization
tn = tsearchn([X,Y,Z], tri, [xy;z]'); % Determine which triangle point is within
IsInside = ~isnan(tn) % Convert to logical vector

% Only retain the points within the cylinder
xy = xy(:, IsInside);
z = z(:, IsInside);
figure(33)
scatter3(xy(1,:), xy(2,:), z)