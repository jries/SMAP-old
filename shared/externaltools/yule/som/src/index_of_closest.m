function i = index_of_closest(x, w)
% I = INDEX_OF_CLOSEST(X, W) takes length-N vector X nad MxN matrix W
% and returns row I of W closest to X.
[ignore, i] = closest(x, w);
