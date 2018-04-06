function w = update_weights(w, u, x_pick, i, t, tmax, ...
                            mu_i, mu_f, sigma_i, sigma_f)
% UPDATE_WEIGHTS  Update weights for Kohonen's Self-Organizing Map.
%     Generally this function should not be called directly, but is called
%     automatically by SOMLEARN.  See SOMLEARN for an explanation of the
%     parameters.

% scale learning paramters by elapsed time                        
tfrac = t / tmax;
mu = scale(mu_i, mu_f, tfrac);
sigma = scale(sigma_i, sigma_f, tfrac);

% udpate the weights, tracking mean weight change
for k = 1:size(w, 1)
    aleph = exp(-sum((u(i,:)-u(k,:)).^2) / (2*sigma^2));
    dw = mu * aleph * (x_pick - w(k,:));
    w(k,:) = w(k,:) + dw; 
end
        