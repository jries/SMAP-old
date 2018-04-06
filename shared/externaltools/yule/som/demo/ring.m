function a = ring(n, ri, ro)
% A = RING(N, RI, RO) generates points of a ring in the unit square of
% N random points,centered at [0.5, 0.5], with inner radius RI and outer 
% radius RO.

[a,c] = disc(n, ro);

c = c(1:size(a, 1),:);

a = a(find(eucdist(a, c) > ri), :);

