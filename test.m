function p = test()
    unitSur = mkSurCom(20, 80, 20, 50);
    unitSur2 = mkSurCom(20, 70, 20, 40);
    figure(1000)
    visualSurCom(unitSur);
    hold on
    visualSurCom(unitSur2);
    hold off
    p=[];
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
    X = reshape(X,[numel(X),1]);
    Y = reshape(Y,[numel(Y),1]);
    Z = reshape(Z,[numel(Z),1]);
    
    Z = Z * unitSur.mainDepth;
    Z = Z + unitSur.root;
end

function callbackfunction(uiobject,data,obj,extradata)
disp('callback')
end

function aftersync_callback(obj)
disp('aftersync')
end