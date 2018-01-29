function [x, y, z] = pointsIntheSpace(depth, dia, thickness, numMol)
    [X, Y, Z] = parabolaCy(depth, dia*2); % make an outer parabola cylinder
    outVol = volume(alphaShape(X,Y,Z)); % calculate the volume of the outer cylinder
    [Xin, Yin, Zin] = parabolaCy(depth-thickness, dia-thickness/2); % make an inner parabola cylinder
    inVol = volume(alphaShape(Xin, Yin, Zin)); % calculate the volume of the inner cylinder
    factor = (outVol-inVol)/(dia*2)^2*depth; % the ratio between the random space and space of interest (space between the two cylinders)
    inputNum = upper(numMol*factor); % calculate how many points need to be sampled
    % numMol = 200
    % Generate random points
    n = ceil(inputNum*1.2); % the total number of random points the 1.2 times of inputNum. The number of points should be higher than numMol.
    xy = rand(2, n).*dia*2-dia; % generate random points at xy
    z = rand(1, n).*depth; % generate random points at z
    
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
    x = x(:,1:numMol)';
    y = y(:,1:numMol)';
    z = z(:,1:numMol)';
    figure(33)
    scatter3(x, y, z)
    figure(34)
    scatter(x, y)
end

function [X, Y, Z] = parabolaCy(depth, dia)
    a = aFinder(depth, dia);
    Y = 0:1:depth;
    X = ((Y-depth)/a).^(1/2);
    [X,Y,Z] = cylinder(X,20);
    X = reshape(X,[numel(X),1]);
    Y = reshape(Y,[numel(Y),1]);
    Z = reshape(Z,[numel(Z),1]);
    Z = Z*depth;
end

function a = aFinder(depth, dia)
    Y = 0:1:depth;
    a = -depth/(dia).^2;
end

