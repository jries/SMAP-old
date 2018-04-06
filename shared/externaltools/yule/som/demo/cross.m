function a = cross(n, t) 
% A = CROSS(N, T) generates points of a cross in the unit square of N random
% points, centered at [0.5 0.5], with thickness T.  

a = rand(n, 2);

lo = 0.5 - t;
hi = 0.5 + t;

a = [block(a, lo, hi, lo, hi); ...
     block(a, hi, 1, lo, hi); ...
     block(a, lo, hi, hi, 1); ...
     block(a, 0, lo, lo, hi); ...
     block(a, lo, hi, 0, lo)];


function a = block(a, xlo, xhi, ylo, yhi)
a = a(find(a(:,1) > xlo & a(:,1) < xhi & a(:,2) > ylo & a(:,2) < yhi), :);

