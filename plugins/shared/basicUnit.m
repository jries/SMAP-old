classdef basicUnit

    properties
        %define class properties if needed
        capDepth
        capShape
        bottomDepth
        bottomShape
        radius
        basicComponent
        
    end
    methods
        function obj=basicUnit(capDepth, capShape, bottomDepth, bottomShape, radius)
            % basicUnit(25,'ellipse', 80, 'bar', 25)
            obj.capDepth = capDepth;
            obj.capShape = capShape;
            obj.radius = radius;
            obj.bottomDepth = bottomDepth;
            obj.bottomShape = bottomShape;
            obj.basicComponent = mkBasicCom(capDepth, capShape, bottomDepth, bottomShape, radius);
        end
    end
end

function  basicComponent = mkBasicCom(capDepth, capShape, bottomDepth, bottomShape, radius)
    %% make a surface component for making a cylinder
    switch capShape
        case 'ellipse'
        a=capDepth; % horizontal radius
        b=radius; % vertical radius
        xO=0; % ellipse centre coordinates
        yO=bottomDepth;
        t=-pi:0.01:pi;
        x=xO+a*cos(t);
        y=yO+b*sin(t);
        idxKept = x<=xO&y>=yO;
        cx = x(idxKept);
        cy = y(idxKept);
    end
   
    if bottomDepth > 0
        switch bottomShape
            case 'bar'
                py = 0:1:bottomDepth;
                px = repelem(-radius, numel(py));
                x = [px cx];
                y = [py cy];
        end
    else
        x = cx;
        y = cy;
    end
    basicComponent = [];
    xy = unique([x;y]','rows','stable');
    yq = 0:(capDepth+bottomDepth);
    xq = interp1(xy(:,2),xy(:,1),yq);
    basicComponent.x = xq;
    basicComponent.y = yq;
end