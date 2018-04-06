function y = mag(x)
% Y = MAG(X) returns magnitude Y = ||X||.

y = sqrt(sum(abs(x).^2, 2));