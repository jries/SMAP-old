function l = pointsIntheSpace(pp)
    [X, Y, Z] = parabolaCy(pp.depth, pp.dia); % make an outer parabola cylinder
    [~, outVol] = convhull(X,Y,Z); % calculate the volume of the outer cylinder
    [Xin, Yin, Zin] = parabolaCy(pp.depth-pp.thickness, pp.dia-pp.thickness/2); % make an inner parabola cylinder
    [~, inVol] = convhull(Xin, Yin, Zin); % calculate the volume of the inner cylinder
    factor = ((pp.dia*2)^2*pp.depth)/(outVol-inVol); % the ratio between the random space and space of interest (space between the two cylinders)
    inputNum = upper(pp.numMol*factor); % calculate how many points need to be sampled
    % numMol = 200
    % Generate random points
    n = ceil(inputNum*1.2); % the total number of random points the 1.2 times of inputNum. The number of points should be higher than numMol.
    xy = rand(2, n).*pp.dia*2-pp.dia; % generate random points at xy
    z = rand(1, n).*pp.depth; % generate random points at z
    
    % Outer cylinder
    % Only retain the points within the cylinder
    tri = delaunayn([X,Y,Z]); % Generate delaunay triangulization
    tn = tsearchn([X,Y,Z], tri, [xy;z]'); % Determine which triangle point is within
    IsInside = ~isnan(tn); % Convert to logical vector
    % Filter the points (in the outer cylinder)
    xy = xy(:, IsInside);
    z = z(:, IsInside);
    
    % Inner cylinder
    % Only retain the points within the cylinder
    tri = delaunayn([Xin,Yin ,Zin]); % Generate delaunay triangulization
    tn = tsearchn([Xin,Yin,Zin], tri, [xy;z]'); % Determine which triangle point is within
    IsInside = ~isnan(tn) % Convert to logical vector
    
    % Filter the points (out of the inner cylinder)
    xy = xy(:, ~IsInside);
    z = z(:, ~IsInside);
    x = xy(1,:);
    y = xy(2,:);
    x = x(:,1:pp.numMol)' + pp.comp;
    y = y(:,1:pp.numMol)' + pp.comp;
    z = z(:,1:pp.numMol)';
    
    switch pp.view
        case 'bottom'
            l.x = x
            l.y = y
        case 'side'
            l.x = x
            l.y = z
        otherwise
    
    end
    
    figure(33)
    scatter3(x, y, z)
    figure(34)
    scatter(x, y)
end

function [X, Y, Z] = parabolaCy(depth, dia)
%% Using a Parabola to create a cylinder
    a = aFinder(depth, dia);
    pY = 0:1:depth;
    pX = ((pY-depth)/a).^(1/2);
    [X,Y,Z] = cylinder(pX, 40); % create a unit cylinder (the range of z is from 0 to 1)
    % get one column vectors of X, Y, and Z
    X = reshape(X,[numel(X),1]);
    Y = reshape(Y,[numel(Y),1]);
    Z = reshape(Z,[numel(Z),1]);
    Z = Z*depth; % assign depth to the cylinder
end

function a = aFinder(depth, dia)
%% This is for the determination of coefficient "a" for a Parabola, using depth (y-intercept) and diameter (x-intercept) as inputs
    Y = 0:1:depth;
    a = -depth/(dia).^2;
end

