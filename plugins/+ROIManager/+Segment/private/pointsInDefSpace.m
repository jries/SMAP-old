function l = pointsInDefSpace(p)
    outUnitSur = mkSurCom(p.outCap, p.outBottom, p.root, p.outDia);
    inUnitSur = mkSurCom(p.inCap, p.inBottom, p.root, p.inDia);
    [Xout, Yout, Zout] = mkCy(outUnitSur);
    [Xin, Yin, Zin] = mkCy(inUnitSur);
    
    totalDepth = p.root+outUnitSur.mainDepth;
    
    [T,R,ZRin] = meshgrid(linspace(0,2*pi,41),linspace(0,p.inBottom/2,16),totalDepth);
    XRin = R.*cos(T);
    YRin = R.*sin(T);
    figure(30)
    scatter3(XRin(:),YRin(:),ZRin(:))
    Xin=[XRin;Xin];
    Yin=[YRin;Yin];
    Zin=[ZRin;Zin];
    
    [T,R,ZRout] = meshgrid(linspace(0,2*pi,41),linspace(0,p.outBottom/2,16),totalDepth);
    XRout = R.*cos(T);
    YRout = R.*sin(T);
    figure(30)
    scatter3(XRout(:),YRout(:),ZRout(:))
    Xout=[XRout;Xout];
    Yout=[YRout;Yout];
    Zout=[ZRout;Zout];

    outSur = scatteredInterpolant(Xout(:), Yout(:), Zout(:));
    outSur.ExtrapolationMethod = 'none';
    outSur.Method = 'natural';
    
    inSur = scatteredInterpolant(Xin(:), Yin(:), Zin(:));
    inSur.Method = 'natural';
    
    
    pointMode = 'full';
    switch pointMode
        case 'random'
        % Generate random points #DEV
        [inK, inVol] = convhull(Xin, Yin, Zin); % calculate the volume of the inner cylinder
        [outK, outVol] = convhull(Xout, Yout, Zout); % calculate the volume of the outer cylinder
        factor = ((p.outDia)^2*totalDepth)/(outVol-inVol); % the ratio between the random space and space of interest (space between the two cylinders)
        inputNum = upper(p.numMol*factor); % calculate how many points need to be sampled
        n = ceil(inputNum*1.1); % the total number of random points the 1.2 times of inputNum. The number of points should be higher than numMol.   % numMol = 200

        xy = rand(2, n).*(p.outDia/2)*2-p.outDia/2; % generate random points at xy
        z = rand(1, n).*totalDepth; % generate random points at z
    
        case 'full'
        % Generate full points
        sizeRange = p.size/2;
        [x y z] = meshgrid((-sizeRange):sizeRange, (-sizeRange):sizeRange, 0:totalDepth);
        xy = [x(:), y(:)]';
        z = z(:)';
    end
    
    % Filter the points
    refZin = inSur(xy(1,:),xy(2,:));
    refZout = outSur(xy(1,:),xy(2,:));
    
    idx = z<refZout&z>refZin&z>=p.root;
    l = [];
    tempL = [x(:), y(:), z'];
    tempL = tempL(idx,:);
    
    tempL(:,1) = tempL(:,1) + sizeRange;
    tempL(:,2) = tempL(:,2) + sizeRange;
    tempL = tempL+1;
      
    outputType = 'image';
    switch outputType
        case 'coordinates'
            %% make the projection by discard infomation of a specified axis
            l = [];
            switch p.viewType.selection
                case 'top'
                    % discard z axis
                    l.x = x(1,1:p.numMol)';
                    l.y = y(1,1:p.numMol)';
                case 'side'
                    % discard y axis, and relabeled oringinal z axis as new y axis
                    l.x = x(1,1:p.numMol)';
                    l.y = z(1,1:p.numMol)';
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
        case 'image'
            %%
            switch p.viewType.selection
                case 'top'
                    % discard z axis
                    proImg = accumarray(int8(tempL(:,1:2)), tempL(:,3),[],@(x) length(x));
                    figure(39)
                    scatter3(tempL(:,1), tempL(:,2), tempL(:,3))
                    rotate3d on
                    figure(31)
                    image(proImg, 'CDataMapping', 'scaled')
                    imwrite(uint8(round(mat2gray(proImg)*255)), [p.folderPath '\' p.imgPath])
                case 'side'
                    % discard y axis, and relabeled oringinal z axis as new y axis

                otherwise
            end
    end

end

function [xy, z] = pointsIn3DS(X, Y, Z, xy, z, filter)
    tri = delaunayn([X,Y,Z]); % Generate delaunay triangulization
    tn = tsearchn([X,Y,Z], tri, [xy;z]'); % Determine which triangle point is within
    switch filter
        case 'inside'
            idx = ~isnan(tn); % Convert to logical vector
        case 'outside'
            idx = isnan(tn); % Convert to logical vector
    end
    
    % Filter the points (in the outer cylinder)
    xy = xy(:, idx);
    z = z(:, idx);
end

function  unitSur = mkSurCom(capDepth, bottomDepth, root, diameter)
    %% make a surface component for making a cylinder
    if capDepth > 0
        % Create a parabola cap
        a = aFinder(capDepth, diameter);
        pY = 0:1:capDepth;
        pX1 = ((pY-capDepth)/a).^(1/2);
        rmBottom = 1;
    else
        pX1=[];
        rmBottom = 0;
    end
    
    if bottomDepth > 0
        % add a bottom column
        pY = 0:1:bottomDepth;
        pX2 = repelem(diameter/2, numel(pY));
    else
        pX2=[];
        rmBottom = 0;
    end
    
    pX = [pX2(1:(numel(pX2)-rmBottom)) pX1];
    
    % set the root
    unitSur = [];
    unitSur.main = pX;
    unitSur.root = root;
    unitSur.mainDepth = capDepth + bottomDepth;
end

function a = aFinder(depth, dia)
%% This is for the determination of coefficient "a" for a Parabola, using depth (y-intercept) and diameter (x-intercept) as inputs
    Y = 0:1:depth;
    a = -depth/(dia/2).^2;
end

function  fig = visualSurCom(unitSur)
    reflection = [-unitSur.main; 0:unitSur.mainDepth];
    surfaceCurve = [reflection(:,end:-1:1), [unitSur.main; 0:unitSur.mainDepth]];
    plot(surfaceCurve(1,:), surfaceCurve(2,:))
end

function  [X, Y, Z] = mkCy(unitSur)
    [X,Y,Z] = cylinder(unitSur.main, 40); % create a unit cylinder (the range of z is from 0 to 1)
       
    Z = Z * unitSur.mainDepth;
    Z = Z + unitSur.root;
end