function [ sigma_x, sigma_y ] = calc_sigma( data, mean_x, mean_y)
    data = squeeze(sum(data, 4));
  %  dipshow(data);
    I = sum(data(:)) / 4;
    size = length(data(:, 1));
    sigma_x = 0;
    sigma_y = 0;
    for ii = 1 : 1 : size
        for jj = 1 : 1 : size
            sigma_x = sigma_x + (data(ii, jj) * (ii - mean_x) ^ 2); 
            sigma_y = sigma_y + (data(ii, jj) * (jj - mean_y) ^ 2);
        end
    end
    sigma_x = sqrt(sigma_x / I);
    sigma_y = sqrt(sigma_y / I);
end

