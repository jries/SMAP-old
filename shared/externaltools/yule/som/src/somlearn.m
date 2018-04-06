function w = somlearn(w, u, x, mu_i, mu_f, sigma_i, sigma_f, t, tmax)
% SOMLEARN  Run one step of learning using Kohonen's Self-Organizing Map.
%    W = SOMLEARN(W, U, X, MU_I, MU_F, SIGMA_I, SIGMA_F, T, TMAX) returns
%    new weights W based on current weights W, output vectors U, and input
%    vectors X.  Remaining parameters specify weight update according to 
%    the formulas:
%
%         mu_t = MU_I + T/TMAX * (MU_F - MU_I)
%    
%         sigma_t = SIGMA_I + T/TMAX * (SIGMA_F - SIGMA_I)
%  
%         aleph(i, k) = e^(-|u_i-u_k|^2 / 2 * sigma_t^2)
%    
%         w_i = w_i + mu_t * aleph(i, k) * (x - w_i) 
%
%     where w_i is the weight vector associate with output vector u_i, and
%     u_k is the output vector that wins the competition on this round.

% randomly choose an input vector x
x_pick = pickrand(x);

% determine the winning output node i closest to x
i = index_of_closest(x_pick, w);

% update weights and track weight changes
w = update_weights(w, u, x_pick, i, t, tmax, mu_i, mu_f, sigma_i, sigma_f);
