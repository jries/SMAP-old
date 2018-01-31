function l = pointsIntheSpace(pp)
    [X, Y, Z] = mkCy(pp.shapeOut, pp.depth, pp.dia); % make an outer parabola cylinder
    figure(30)
    scatter3(X,Y,Z)
    rotate3d on
    [~, outVol] = convhull(X,Y,Z); % calculate the volume of the outer cylinder
    
    % Determine the thickness of the top
    if regexp(pp.shapeIn, 'NT$')>1
        zThickness = 0;
    else
        zThickness = pp.thickness;
    end
    [Xin, Yin, Zin] = mkCy(pp.shapeIn, pp.depth-zThickness, pp.dia-pp.thickness); % make an inner parabola cylinder
    figure(31)
    scatter3(Xin,Yin,Zin)
    rotate3d on
    [~, inVol] = convhull(Xin, Yin, Zin); % calculate the volume of the inner cylinder
    factor = ((pp.dia*2)^2*pp.depth)/(outVol-inVol); % the ratio between the random space and space of interest (space between the two cylinders)
    inputNum = upper(pp.numMol*factor); % calculate how many points need to be sampled
    % numMol = 200
    % Generate random points
    n = ceil(inputNum*1.1); % the total number of random points the 1.2 times of inputNum. The number of points should be higher than numMol.
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
    x = x(:,1:pp.numMol)' + pp.comp; % pp.comp ensures the sturctural sizes of two channels are the same
    y = y(:,1:pp.numMol)' + pp.comp; % the same
    z = z(:,1:pp.numMol)';
    
    %% make the projection by discard infomation of a specified axis
    switch pp.view
        case 'bottom'
            % discard z axis
            l.x = x
            l.y = y
        case 'side'
            % discard y axis, and relabeled oringinal z axis as new y axis
            l.x = x
            l.y = z
        otherwise
    
    end
    
    % visualization of the structures
    figure(33)
    scatter3(x, y, z)
    rotate3d on
    figure(34)
    scatter(x, y)
    figure(35)
    scatter(x, z)
end

function [X, Y, Z] = mkCy(shape, depth, dia)
switch shape
    case 'parabola'
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
    case 'columnRT'
    %% create a round-top cylinder
        depthTop = 20; % the range of curved region (from top)
        depthBottom = depth - depthTop; % the range of column region (from bottom)
        
        % top region
        a = aFinder(depthTop, dia);
        pY = 0:1:depthTop;
        pX = ((pY-depthTop)/a).^(1/2);
        [Xtop,Ytop,Ztop] = cylinder(pX, 40); % create a unit cylinder (the range of z is from 0 to 1)
        % get one column vectors of X, Y, and Z
        Xtop = reshape(Xtop,[numel(Xtop),1]);
        Ytop = reshape(Ytop,[numel(Ytop),1]);
        Ztop = reshape(Ztop,[numel(Ztop),1])*depthTop+depthBottom;
        
        % bottom region
        pY = 0:1:depthBottom;
        pX = repelem(dia, numel(pY));
        [Xbottom,Ybottom,Zbottom] = cylinder(pX, 40);
        Xbottom = reshape(Xbottom,[numel(Xbottom),1]);
        Ybottom = reshape(Ybottom,[numel(Ybottom),1]);
        Zbottom = reshape(Zbottom,[numel(Zbottom),1]);
        Zbottom = Zbottom*depthBottom; % assign depth to the cylinder
        
        X = [Xbottom' Xtop']';
        Y = [Ybottom' Ytop']';
        Z = [Zbottom' Ztop']';
    case 'columnNT'
    %% create a cylinder with no top
        pY = 0:1:depth;
        pX = repelem(dia, numel(pY));
        [X,Y,Z] = cylinder(pX, 40);
        X = reshape(X,[numel(X),1]);
        Y = reshape(Y,[numel(Y),1]);
        Z = reshape(Z,[numel(Z),1]);
        Z = Z*depth; % assign depth to the cylinder

    otherwise
end

end

function a = aFinder(depth, dia)
%% This is for the determination of coefficient "a" for a Parabola, using depth (y-intercept) and diameter (x-intercept) as inputs
    Y = 0:1:depth;
    a = -depth/(dia).^2;
end

