function G = wanderlust(data, Options)
%    G = wanderlust( data, Options) 
%
%    Options                    : structure for specifying algorithm options which can contain
%                                   zero or more of the following fields:
%
%      [k]                      : truncate the neighbors to random k out of l 
%      [l]                      : size of neighborhood l>k closest points
%      [s]                      : index to the starting point in the data
%      [num_graphs]             : number of repeats - random selections of
%                                   k out of l.
%      [num_landmarks]          : number of waypoints\landmarks <- TODO
%                                 ensure logic if someone specifies predefined landmarks
%      [verbose]                : messages are printed to Matlab's sto
%      [metric]                 : string - distance metric for constructing 
%                                   the nearest neighbor graphs
%      [voting_scheme]          : CYCLER: How to weigh each point's contribution
%           'uniform'           -  
%           'exponential'       - 
%           'linear'            - 
%           'inverse-linear'    - 
%           TOOD return k_power -
%      [band_sample]            : CYCLER: true\false landmarks are subsamples at equidistanced
%                                   bands away from the start point. using 
%                                   the points shortest path
%                                   distance over the graph.
%      [partial_order]          : an array of indices in data. the indices 
%                                point to landmarks and their order is forced 
%                                in the final output. 
%      [flock_landmarks]        : CYCLER: number of how many times to use a median
%      [plot_data]              : 2D data matrix to plot results to.
%                                  filter on landmarks
%      [snn]                    : shared nearest neighbor
%      [ann]                    : TODO - finalize - using Approx NN
%      [search_connected_components] :TODO add the option to cancel
%      
%
% run the wonderlust algorithm on data. generate num_graphs klNN graphs; starting point is s, choose num_landmarks
% random landmarks. weigh landmark-node pairs by distance to the power of power_k.
%
% (use distance as distance metric)
%
% (alternatively data can be lnn graph)
% 

% set up return structure
G.landmarks = [];
G.T = []; % traj
G.B = []; % branch

% Algorithm defaults
G.Opts.metric = 'euclidean';
G.Opts.k = 20;
G.Opts.l = 30;
G.Opts.num_graphs = 10;
G.Opts.s = randi(30, 1);% G.Opts.s = randi(size(data, 1), 1);
G.Opts.num_landmarks = 30; %G.Opts.num_landmarks = min(size(data, 1), 100);
G.Opts.verbose = true;
G.Opts.branch = false;
G.Opts.partial_order = [];
G.Opts.voting_scheme = 'exponential';
G.Opts.band_sample = false;
G.Opts.flock_landmarks = 2;
G.Opts.search_connected_components = true;
G.Opts.plot_landmark_paths = false;
G.Opts.plot_data = [data(:,1) data(:,2)];

fn = fieldnames(Options);
%while 0
for j=1:length(fn)
    name = fn{j};
    value = getfield(Options, name);
    
    if     strcmpi(name,'metric')           G.Opts.metric = value;
    elseif strcmpi(name,'k')                G.Opts.k = value;
    elseif strcmpi(name,'l')                G.Opts.l = value;
    elseif strcmpi(name,'num_graphs')       G.Opts.num_graphs = value;
    elseif strcmpi(name,'s')                G.Opts.s = value;
    elseif strcmpi(name,'num_landmarks')    G.Opts.num_landmarks = value;
    elseif strcmpi(name,'verbose')          G.Opts.verbose = value;
    elseif strcmpi(name,'branch')           G.Opts.branch = value;
    elseif strcmpi(name,'partial_order')   	G.Opts.partial_order = value;
    elseif strcmpi(name,'band_sample')      G.Opts.band_sample = value; 
    elseif strcmpi(name,'flock_landmarks')  G.Opts.flock_landmarks = value; 
    elseif strcmpi(name,'plot_data')        G.Opts.plot_data = value; 
    elseif strcmpi(name,'plot_landmark_paths') G.Opts.plot_landmark_paths = value; 
%     else   fprintf('Cycler: invalid option "%s" ignored.\n', name);
    end
end

% build lNN graph
if issparse( data ) 
    if G.Opts.verbose 
        disp 'using prebuilt lNN graph';
    end
    lnn = data;
else
    % after normalization ndata may have duplicates
    [u,I,J] = unique(data, 'rows', 'first');
    G.hasDuplicates = size(u,1) < size(data,1);

    if G.hasDuplicates
        disp 'warn: data has duplicate values';
        
        G.J = J;
        G.Opts.s  = find(ismember(u, data(G.Opts.s, :), 'rows'));
        data = u;
    end

    if G.Opts.verbose 
        disp 'building lNN graph';
    end

    tic
	lnn = parfor_spdists_knngraph( data, G.Opts.l,...
        'distance', G.Opts.metric,...
        'chunk_size', min(500, size( data, 1 )),... % TODO: parameterize and add opt for ppl without PC toolbox
        'verbose', G.Opts.verbose );
    
    sprintf('lnn computed: %gs', toc/60);
end

% generate klNN graphs and iteratively refine a trajectory in each
for graph_iter = 1:G.Opts.num_graphs
	if( G.Opts.verbose )
		fprintf( 1, 'iter #%d:\n', graph_iter );
	end

	iter_t = tic;

	% randomly generate a klNN graph
	if( G.Opts.verbose )
		fprintf( 1, 'entering knn graph: ' );
	end

	klnn = spdists_klnn( lnn, G.Opts.k, G.Opts.verbose );
 	klnn = spdists_undirected( klnn ); % TODO consider removing 
%     i = graph(klnn);
%     h = plot(i);
%     highlight(h,1101:1380, 'Nodecolor', 'r')
%     figure(10024)
%     h
    
	if( G.Opts.verbose )
		fprintf( 1, ' done (%3.2fs)\n', toc( iter_t ) );
		fprintf( 1, 'entering trajectory landmarks: ' );
	end

	% run traj. landmarks
	[ traj, dist, iter_l, RNK,paths_l2l ] = trajectory_landmarks( klnn,data, G);
    G.landmarks(graph_iter, :) = iter_l;
    G.traj = traj;
    G.dist = dist;
    
    if G.Opts.verbose
        fprintf( 1, ' done (%3.2fs)...\n', toc( iter_t ) );
    end

	% calculate weighed trajectory    
    if strcmpi(G.Opts.voting_scheme, 'uniform')
        W(:, :) = 1;
    elseif strcmpi(G.Opts.voting_scheme, 'exponential')
        sdv = mean ( std ( dist) )*3;
        W = exp( -.5 * (dist / sdv).^2);
    elseif strcmpi(G.Opts.voting_scheme, 'linear')
        W = repmat(max(dist), size( dist, 1 ), 1) - dist;
    end
    
    % The weghing matrix must be a stochastoc operator
    W = W ./ repmat( sum( W ), size( W, 1 ), 1 );
        
    t( 1,:)  = traj(1,:);
	t( 2, : ) = sum( traj .* W );
    
    if any(isnan(t(2,:)))
        fprintf('error: there are nans in trajectory (unreachable points) \n');
        return;
    end
    
	% iteratively realign trajectory (because landmarks moved)
	converged = 0; user_break = 0; realign_iter = 2;
	while( ~converged && ~user_break)
		realign_iter = realign_iter + 1;

		traj = dist;
		for idx = 1:size( dist, 1 )
			% find position of landmark in previous iteration
			idx_val = t( realign_iter - 1, iter_l( idx ) );
			% convert all cells before starting point to the negative
			before_indices = find( t( realign_iter - 1, : ) < idx_val );
			traj( idx, before_indices ) = -dist( idx, before_indices );
			% set zero to position of starting point
			traj( idx, : ) = traj( idx, : ) + idx_val;
        end
        
		% calculate weighed trajectory
		t( realign_iter, : ) = sum( traj .* W );

		% check for convergence
        corr( t( realign_iter, : )', t( realign_iter - 1, : )' )
		converged = corr( t( realign_iter, : )', t( realign_iter - 1, : )' ) > 0.9995;
        
        if (mod(realign_iter,500)==0)
            user_break = true;
        end
	end
	fprintf( 1, '%d realignment iterations, ', realign_iter );

	% save final trajectory for this graph
    if G.hasDuplicates        
        G.T(graph_iter, :) = t(realign_iter, G.J);
        G.s  = find(ismember(data(G.J,:), data(G.Opts.s, :), 'rows'));
    else
        G.T(graph_iter, :) = t(realign_iter, :);
    end
    
    if ( G.Opts.verbose )
		toc( iter_t );

		fprintf( 1, '\n' );
	end
end
end


function spdists = spdists_klnn( spdists, k, verbose )
	% spdists = spdists_klnn( spdists, k, verbose )
	%
	% given a lNN graph spdists, choose k neighbors randomly out of l for each node

	remove_edges = [];

	for idx = 1:length( spdists )
		% remove l-k neighbors at random
		neighs = find( spdists( :, idx ) );
		l = length( neighs ); % count number of neighbors
		remove_indices = neighs( randsample( length( neighs ), l - k ) );
		idx_remove_edges = sub2ind( size( spdists ), remove_indices, ones( l - k, 1 ) * idx );
		remove_edges = [ remove_edges; idx_remove_edges ];

		if( verbose )
			if( mod( idx, 50000 ) == 0 )
				fprintf( 1, '%3.2f%%', idx / length( spdists ) * 100 );
			elseif( mod( idx, 10000 ) == 0 )
				fprintf( 1, '.' );
			end
		end
	end

	spdists( remove_edges ) = 0;
    end

	function [ traj, dist, l, RNK,paths_l2l ] = trajectory_landmarks( spdists,data, G)
	% [ traj, dist, l ] = trajectory_landmarks( spdists, s, n, verbose )
	%
	% calculate the trajectory score of each point in spdists.
	%
	% s: list of indices of possible starting points. one of these points will be used to generate a reference
	% trajectory; the landmark shortest paths will be aligned to this reference.
	% n: list of landmark indices to use; or, alternatively, the number of landmarks to choose randomly from all
	% points.
	%
	% traj is a |n|x|spdists| matrix, where row i is the aligned shortest path from landmark i to each other point.
	% dist is a |n|x|spdists| matrix, where row i is the shortest path from landmark i to each other point. l is
	% the list of landmarks, l(1) is the starting point.

    RNK = zeros(size(data, 1), 1);
    n = G.Opts.num_landmarks;

    if( length( G.Opts.s ) > 1 )
		% if given a list of possible starting points, choose one. TODO move
		% to beginning of algorithm!!!!
		G.Opts.s = randsample( G.Opts.s, 1 );
	end

	if( length( n ) == 1 )
        [dists, paths, ~] = graphshortestpath( spdists, G.Opts.s,'METHOD','Dijkstra');%, 'directed', true);
        
        % if not given landmarks list, decide on random landmarks        
        if (G.Opts.band_sample)
            n_opts = [];
            num_jumps_arr = cellfun(@(x)numel(x), paths);
            max_jumps = max(num_jumps_arr);
            max_dist = max(dists);
            for (prc = .95:-.10:.15)
                band = find(dists>=(prc-.1)*max_dist & dists <=prc*max_dist); % & num_jumps_arr >= floor((prc-.05)*max_jumps) & num_jumps_arr <= ceil(prc*max_jumps));
                if length(band)> (n- 1 - length(G.Opts.partial_order))
                    n_opts = [n_opts randsample( band, n - 1 - length(G.Opts.partial_order), true )];
                else
                    n_opts = [n_opts band];
                end
            end
            n = randsample( n_opts, n - 1 - length(G.Opts.partial_order) );
        else
            n = randsample( 1:size(data,1), n - 1 - length(G.Opts.partial_order) );
        end
        
        % flock landmarks 
        if (G.Opts.flock_landmarks > 0)
        for k=1:G.Opts.flock_landmarks
            [IDX, ~] = knnsearch(data, data(n, :), 'distance', G.Opts.metric, 'K', 20);     
            for i=1:numel(n)
                n(i) = knnsearch(data, median(data(IDX(i, :), :)), 'distance', G.Opts.metric); 
            end
        end
        end
    end

    partial_order = [G.Opts.s;G.Opts.partial_order(:)]; % partial_order includes start point
	l = [ partial_order; n(:) ]; % add extra landmarks if user specified
    
    % prune out weak edges
    prune = false;
    if (prune)
        res = 256;
        [band, density, x, y] = kde2d(data, res);

            % maybe remove edges jumping over low density regions?
            % BW1 = edge(density,'canny', .001);
            % imshow(BW1)
            % set(gca,'YDir','normal')

            %remove bad edges if the path to a landmark has a suspicious jump
        mapX = round( interp1(diag(x), 1:res, data(:, 1), 'linear', 'extrap') );
        mapY = round( interp1(diag(y), 1:res, data(:, 2), 'linear', 'extrap') );
        recalc = true;
        cou = 0;
        while recalc && cou < 30
            cou = cou+1;
            [dist( 1, : ), paths, ~] = graphshortestpath( spdists, s);%, 'directed', false );
            recalc = false;
            for pathidx=2:numel(l) %iterate over path to each landmark
                l_path = paths{l(pathidx)};
                path_jumpsX = abs(mapX(l_path(2:end))-mapX(l_path(1:end-1)));
                path_jumpsY = abs(mapY(l_path(2:end))-mapY(l_path(1:end-1)));
                path_jumps = path_jumpsX + path_jumpsY;
                bad_nodes = find(path_jumps > (mean(path_jumps) + 4*std(path_jumps)));
                if any(bad_nodes)
                    disp(sprintf('removing %g bad connections\n', numel(bad_nodes)));
                    spdists(sub2ind(size(spdists), l_path(bad_nodes), l_path(bad_nodes+1))) = 0;
                    spdists(sub2ind(size(spdists), l_path(bad_nodes+1), l_path(bad_nodes))) = 0;
                    recalc = true;
                end
            end
        end
    end    

	% calculate all shortest paths
    paths_l2l = cell(length(l));
    for li = 1:length( l )
        [dist( li, : ), paths, ~] = graphshortestpath( spdists, l( li ),'METHOD','Dijkstra');%, 'directed', false );
        paths_l2l(li) = {paths(l)};
        unreachable = (dist(li,:)==inf);
        while (any(unreachable) && G.Opts.search_connected_components)
            fprintf(['\n Warning: %g were unreachable. try increasing l'...
                'or k.Your data is possibly non continous, ie '...
                'has a completely separate cluster of points.'...
                'Wanderlust will roughly estimate their distance for now \n'],...
                sum(unreachable));
            % find closest unreachable point to reachable points.
            % connect it on the spdists. continue iteratively.
            unreachablei = find(unreachable);
            reachablei = find(~unreachable);
            while ~isempty(unreachablei)
                
                [idx, d] = knnsearch(data(unreachablei, :), data(reachablei, :));
                closest_reachable = d==min(d);
                
                %add connection to spdists
                spdists(reachablei(closest_reachable),...
                    unreachablei(idx(closest_reachable))) = min(d);
                spdists(unreachablei(idx(closest_reachable)),...
                    reachablei(closest_reachable)) = min(d);
                % move points from unreachable list to reachable
                reachablei(end+1:end+length(find(closest_reachable))) = ...
                    unreachablei(idx(closest_reachable));
                unreachablei(idx(closest_reachable)) = [];
                
                if sum(unreachable) > 100
                    break;
                end
            end
            [dist( li, : ), paths, ~] = graphshortestpath( spdists, l( li ),'METHOD','Dijkstra');%, 'directed', false );
            paths_l2l(li) = {paths(l)};
            unreachable = (dist(li,:)==inf);
        end
        
        if( G.Opts.verbose )
            fprintf( 1, '.' );
        end
    end
    
    
    % adjust paths according to partial order by redirecting
    nPartialOrder = length(partial_order);
    for radius = 1:nPartialOrder 
        for landmark_row = 1:nPartialOrder
            if (landmark_row + radius <= nPartialOrder)
                a = landmark_row;
                b = landmark_row + (radius-1);
                c = landmark_row + radius;
                dist(a, partial_order(c)) = dist(a, partial_order(b)) + dist(b, partial_order(c));
            end
            if (landmark_row - radius >= 1)
                a = landmark_row;
                b = landmark_row - (radius-1);
                c = landmark_row - radius;
                dist(a, partial_order(c)) = dist(a, partial_order(b)) + dist(b, partial_order(c));
            end
        end
    end

	% align to dist_1 - this for loop refers to partial order stuff
	traj = dist;
    for idx = 2:length(partial_order)
        [~, closest_landmark_row] = min(dist); %closest landmark will determine directionality
        traj(idx, closest_landmark_row < idx) = -dist(idx, closest_landmark_row < idx);
        traj( idx, : ) = traj( idx, : ) + dist( 1, l( idx ) );
    end
    
    % This is the actual align for regular wanderlust
    if length( l ) > length(partial_order)
        for idx = length(partial_order)+1:length( l )
            % find position of landmark in dist_1
            idx_val = dist( 1, l( idx ) );
            % convert all cells before starting point to the negative
            before_indices = find( dist( 1, : ) < idx_val );
            traj( idx, before_indices ) = -dist( idx, before_indices );
            % set zero to position of starting point
            traj( idx, : ) = traj( idx, : ) + idx_val;
        end
    end
    if (G.Opts.plot_landmark_paths)
        plot_landmark_paths(G.Opts.plot_data, paths_l2l, l);
    end
end

function plot_landmark_paths(data, paths, l)
    figure('Color',[1 1 1]);
    scatter (data(:, 1), data(:, 2));
    hold on;
    plot(data(l(1), 1), data(l(1), 2), 'Xr');
    scatter(data(l, 1), data(l, 2), 20*ones(numel(l), 1), 'or');
    for p=1:numel(paths)
        pathsp = paths{p};
        for q=1:numel(pathsp)
            pathq=pathsp{q};
            plot(data(pathq, 1), data(pathq, 2), 'k-');
        end
    end
end