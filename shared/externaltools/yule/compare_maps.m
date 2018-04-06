function cost= compare_maps(map1,map2)
% This function is using t-SNE's cumputing of matrix P and matrix Q in a 
% low-dimensional maps.

    P= compute_Pvalues(map1);
    Q= compute_Pvalues(map2);
    
    const = sum(P(:) .* log(P(:)));                                        % constant in KL divergence
    cost = const - sum(P(:) .* log(Q(:)));                                 %Kullback-Leibler divergence of the maps
end

function value = compute_Pvalues(map) 
% Computing joint probability that point i and j are neighbors
    n = size(map, 1);                                                      % Initialize number of instances
    sum_map = sum(map .^ 2, 2);                                            % make sure the value sum to one
    num = 1 ./ (1 + bsxfun(@plus, sum_map, ...
        bsxfun(@plus, sum_map', -2 * (map * map'))));                      % Student-t distribution
    num(1:n+1:end) = 0;                                                    % set diagonal to zero
    value = max(num ./ sum(num(:)), realmin); 
end


