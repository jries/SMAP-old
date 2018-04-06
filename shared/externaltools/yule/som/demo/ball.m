function [a,c] = ball(n, r)
% A = BALL(N, R) generates points of a ball in the unit cube of
% N random points,centered at [0.5, 0.5, 0.5], with radius R.

a = rand(n, 3);

c = repmat([0.5 0.5 0.5], n, 1);

a = a(find(eucdist(a, c) < r), :);

