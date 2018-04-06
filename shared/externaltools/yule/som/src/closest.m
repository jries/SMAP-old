function [y, i] = closest(x, ys)
% [Y, I] = CLOSEST(X, YS) returns point Y and index I of closest point in
% YS to X.

xs = repmat(x, size(ys, 1), 1);

d = eucdist(xs, ys);

s = sum(d, 2);

[ignore, i] = min(s);

y = ys(i,:);