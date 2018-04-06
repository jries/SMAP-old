function [a,w]  = som(x, n, tmax, seed)
% SOM Kohonen's two-dimensional Self-Organizing Map
%     A  = SOM(X, N, TMAX) returns an N-by-N matrix A of output nodes 
%     learned as an unsupervised map from input vectors X, using Kohonen's 
%     Self-Oganizing Map algorithm.  Matrix A contains indices of X. TMAX 
%     is the number of time steps to to run.   Empirically-determined 
%     learning  parameters are hard-coded in the source code and can be 
%     modified.
%
%     [A,W]  = SOM(X, N, TMAX) also returns the N^2 weights on A.
%
%     SOM(X, N, TMAX, SEED) supports seeding the random-number 
%     generator for reproducible results.
%
%     See also SOMLEARN, PLOT_NODES.

% learning parameters

MU_I        = 0.5;          % learning rate: initial
MU_F        = 0.1;          %              : final

SIGMA_I     =  3.0e0;       % attraction between points : initial
SIGMA_F     =  1.0e-1;      %                           : final

% seed random-number generator if specified
if nargin > 7
    rand('state', seed)
end

% set up NxN grid
[o(:,1), o(:,2)] = ind2sub([n n], 1:n^2);

% create random initial weights
w = rand(n^2, size(x,2));

% run SOM learning for specified number of steps
for t = 1:tmax
    w = somlearn(w, o, x, MU_I, MU_F, SIGMA_I, SIGMA_F, t, tmax);
end

% return grid of input indices
a = zeros(n,n);
for i = 1:size(o,1)
    j = o(i,1);
    k = o(i,2);
    a(j,k) = index_of_closest(w(i,:), x);
end



