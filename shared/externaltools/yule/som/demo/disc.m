function [a,c] = disc(n, r)
% A = DISC(N, R) generates points of a disc in the unit square of
% N random points,centered at [0.5, 0.5], with radius R.

a = rand(n, 2);

c = repmat([0.5 0.5], n, 1);

a = a(find(eucdist(a, c) < r), :);

