classdef LocTransformN<interfaces.LocTransformN0
    %extends pure new transform by legacy functions from LocTransform
    properties
        tinfo
    end
    
    methods
        function [xo,yo,zo]=transformCoordinatesFwd(obj,x,y,z) 
            ci=horzcat(x,y);
            if nargin>3
                ci=horzcat(ci,z);
            end
            co=obj.transformToTarget(2,ci,'nm');
            xo=co(:,1);yo=co(:,2);
            if nargin>3
                zo=co(:,3);
            end
        end
        function [xo,yo,zo]=transformCoordinatesInv(obj,x,y,z) 
            ci=horzcat(x,y);
            if nargin>3
                ci=horzcat(ci,z);
            end
            co=obj.transformToReference(2,ci,'nm');
            xo=co(:,1);yo=co(:,2);
            if nargin>3
                zo=co(:,3);
            end
        end
        function ind=getRef(obj,x,y)
            ci=horzcat(x,y);
            ind=obj.getPart(1,ci,'nm');
        end
        function tinfo=get.tinfo(obj)
            tinfo.mirror.targetmirror=obj.mirror;
            tinfo.separator=[obj.info{1}.xrange(end)*obj.info{1}.cam_pixnm(1) obj.info{1}.yrange(end)*obj.info{1}.cam_pixnm(end)];
        end
    end
end



