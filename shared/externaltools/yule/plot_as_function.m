function plot_as_function(x, Y, varargin)
% plots Y (NxM) as a function of vector x (N length).  
% an average of points in Y is computed for 'num_locs' along vector x.
%
% 'labels': a string cell array to enter for the legend
%
% specify 'avg_type': default 'gaussian_var'
%
% 'sliding': a fixed averaged is computed within a sliding window. takes
% into account points within two windows to each side. so let's say 100
% locations would mean each location is averaging points from two windows to
% the right and to the left
%
% 'linear': same as sliding the each point's contibution dimishes linearly
% 
% 'squared': same as sliding the each point's contibution dimishes
% polynomialy
%
% 'gaussian': all points contibute to the weighted average but beyond 2
% standard deviations the contibution is insignificant.
%
% 'gaussian_var': same as gaussian but the standard deviation is variable
% as a function of the normalized density of the points along x. so in
% saprser locations on x, we use a larger standard deviation to let farther
% points contribute.
%
%
% 'show_error': true or false, shows error bars along the curve, default
% 'false'
% 
% Michelle Tadmor, Columbia University, 2013

clear persistent_ksdensity;
legend_flag = false;
avg_type = 'gaussian_var';
num_locs = 100;
show_error = false;
normalize = true;
rank = false;
make_distribution_movie = false;
smoothness_factor = 0.5;
svGolay = true;
control_density = false;
matPatchColors = [0.75 0.75 0.75; 0.6 0.6 0.6; 0 0 1];
branch = zeros(1, numel(x));

for i=1:length(varargin)-1
    if(strcmp(varargin{i},'num_locs'))
        num_locs = varargin{i+1};
    elseif(strcmp(varargin{i},'avg_type'))
        avg_type = varargin{i+1};
    elseif(strcmp(varargin{i},'labels'))
        legend_flag = true;
        labels = varargin{i+1};
    elseif(strcmp(varargin{i},'show_error'))
        show_error = varargin{i+1};  
    elseif(strcmp(varargin{i},'normalize'))
        normalize = varargin{i+1};
    elseif(strcmp(varargin{i},'rank'))
        rank = varargin{i+1};
        if rank
            control_density = false;
        end
    elseif(strcmp(varargin{i},'svGolay'))
        svGolay = varargin{i+1};
    elseif(strcmp(varargin{i},'smooth'))
        smoothness_factor = varargin{i+1};
    elseif(strcmp(varargin{i},'branch'))
        branch = varargin{i+1};
    elseif(strcmp(varargin{i},'movie'))
        make_distribution_movie = true;
        movie_name = varargin{i+1};
    end
end

if (rank)
    x = tiedrank(x);
end    

weights = zeros(num_locs, length(x));

% compute a weight for each value (data point), at each plot location
for i=1:num_locs
	weights(i, :) = compute_weights(x, (i/num_locs)*range(x)+min(x), avg_type, smoothness_factor);
end

real_weights = weights;
for bri=unique(branch)'

% Compute weighted averages at each location
X = linspace(min(x), max(x), num_locs);
Y_vals = weights*Y./repmat(sum(weights, 2), 1, size(Y, 2));

if control_density
    dens = sum(weights, 2)';
    dens = dens./max(dens);
    dens = floor(dens*10);
    X = linspace(min(x), max(x), num_locs+sum(dens));
    Yrow = 0;
	Y_vals_stretched = zeros(0, size(Y_vals, 2));
    for densi=1:numel(dens)
        Yrow = Yrow+1;
        copies = 0;
        while copies<=dens(densi)
            Y_vals_stretched(end+1, :) = Y_vals(Yrow, :);
            copies = copies+1;
        end
    end
    Y_vals = Y_vals_stretched;
	
    % smooth the streching
    for col=1:size(Y_vals, 2)
        Y_vals(:, col) = smooth(X, Y_vals(:, col),sqrt(num_locs*2), 'sgolay');
    end          

end

if (svGolay && ~show_error) % do not smooth when showing variance
    try 
        for col=1:size(Y_vals, 2)
            Y_vals(:, col) = smooth(X, Y_vals(:, col),sqrt(num_locs*2), 'sgolay');
        end 
    catch 
        disp 'warn: savitzky-golay smoothing is not supported - perhaps you are missing the image smoothing matlab package - no big deal';
    end
end

% Y_vals(2:end,:) = diff(Y_vals);
% Y_vals(1,:) = 0;
Y_vals_raw = Y_vals;
if (normalize)
    % we want to normalize to [0 1]
    mins = prctile(Y_vals, 0, 1);
    Y_vals = bsxfun(@minus, Y_vals, mins);

    rngs = prctile(Y_vals, 100, 1);
    Y_vals = bsxfun(@rdivide, Y_vals, rngs);
end

matColors = distinguishable_colors(size(Y, 2));
set(gca, 'ColorOrder', matColors);
set(0, 'DefaultAxesColorOrder', matColors);

plot(X, Y_vals(:, 1), 'Color', matColors(1, :));    
if (size(Y, 2)> 1)
    for col=2:size(Y, 2)
        hold on;
        plot(X, Y_vals(:, col), 'Color', matColors(col, :));        
    end
end
if (show_error)

    % compute variace along X (symmetically)
    for i=1:num_locs       
        %symmetrical
        Y_errs = bsxfun(@minus,Y,(Y_vals_raw(i, :)));
        
        M = sum(weights(i, :)~=0);
        s = (M-1)/M;
        w_sum = sum(weights(i, :));
        
        Y_valerrs(i, :) = sqrt((weights(i, :)*((Y_errs).^2))/(s*w_sum));
    end
    
	if (normalize)
        Y_valerrs = bsxfun(@rdivide,Y_valerrs,rngs);
	end
     
    % plot the variance as a pretty translucent cloud around line
    for yi=1:size(Y, 2)
        
        % plot light grey background first for variance
        fill( [X, fliplr(X)], [(Y_vals(:, yi)-Y_valerrs(:, yi)./2)', fliplr((Y_vals(:, yi)+Y_valerrs(:, yi)./2)')],...
            matColors(yi,:),'linestyle','none','facealpha',.5);
    end
end

end
if ~(rank || control_density)    
    try
    dens = sum(weights, 2)';
    ca = axis;
    hold on;
    imagesc(X, ca(3):0.04:ca(3)+.05, dens, [0, max(dens)]);
    colorbar;
    axis(ca);
    catch e
    disp(getReport(e,'extended'));
    end
end

if (legend_flag)
    legend(remove_repeating_strings(labels), 'Interpreter', 'none');
end

end

function weights = compute_weights(points, loc, type, factor)
    
    range = quantile(points, .98) - quantile(points, .02);
    min_std_dev = factor*.18*range; % minimum std_dev for dense regions
    max_std_dev = .19*range; % max std_dev for sparse regions  
    linear_slope = 10/range;
    
    if strcmpi('sliding', type) %set '1's on the indices in the windows 
        weights = (points < (loc + 2*min_std_dev)) & ...
                  (points > (loc - 2*min_std_dev));
    
    elseif strcmpi('linear', type)
        weights = 1 - linear_slope*(abs(points - loc));
        weights(weights<0) = 0;
        weights = weights/sum(weights);       
    
    elseif strcmpi('squared', type)
        weights = 1 - ((linear_slope*(points - loc)).^2);
        weights(weights<0) = 0;
        weights = weights/sum(weights);
             
    elseif strcmpi('gaussian_var', type)
        [f, xi] = persistent_ksdensity(points);
    	d = f(minind(abs(xi-loc)))/max(f);

        std_dev = d*(min_std_dev)+(1-d)*max_std_dev;

        weights = ((2*pi*(std_dev^2))^(-1))*exp(-.5*((points - loc)/std_dev).^2);
    
    else % default is strcmpi('gaussian', type)
        weights = ((2*pi*(min_std_dev)^2)^(-1))*exp(-.5*((points - loc)/min_std_dev).^2);
    end
% 	weights = mynormalize(weights(:), 100);
end

function ind = minind(x)
    [~, ind] = min(x);
end


function [f, xi] = persistent_ksdensity(points)

    persistent f_p;
    persistent xi_p;    

    if isempty(f_p)
        [f_p, xi_p] = ksdensity(points);
    end
    f = f_p;
    xi = xi_p;
end