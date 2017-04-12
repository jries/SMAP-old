classdef LocTransform<handle
    properties
        transform %forward .inversie
        tinfo
        transformz;
    end
    methods

        function findTransform(obj,xreference,yreference,xtarget,ytarget,transform)
            %XXX make sure reference and target are not exchanged 
            %matrix becomes singular for too large coordinates. Use um
            %instead of nm
            xreference=xreference/1000;yreference=yreference/1000;
            xtarget=xtarget/1000;ytarget=ytarget/1000;
            switch transform.type
                case {'lwm','polynomial'}
                    obj.transform.inverse = fitgeotrans(double([xtarget ytarget]), [double(xreference) double(yreference)],transform.type,transform.parameter);
                    obj.transform.forward = fitgeotrans( [double(xreference) double(yreference)],double([xtarget ytarget]),transform.type,transform.parameter);
                otherwise             
                     obj.transform.inverse= fitgeotrans(double([xtarget ytarget]),[double(xreference) double(yreference)],transform.type);
                     obj.transform.forward= fitgeotrans([double(xreference) double(yreference)],double([xtarget ytarget]),transform.type);
            end
        end
        
        function findTransformZ(obj,xreference,yreference,zreference,xtargeti,ytargeti,ztarget,transform)
            %do 2D transform first;
%             obj.findTransform(xreference,yreference,xtargeti,ytargeti,transform);
%             transform x,y target
            [xtarget,ytarget]=transformCoordinatesFwd(obj,xtargeti,ytargeti);

            %XXX make sure reference and target are not exchanged 
            %matrix becomes singular for too large coordinates. Use um
            %instead of nm
            xreference=xreference/1000;yreference=yreference/1000;
            xtarget=xtarget/1000;ytarget=ytarget/1000;
            zreference=zreference/1000;
            ztarget=ztarget/1000;
            
            %only affine3d possible      
            obj.transformz.inverse= findAffineTransformZ(double([xtarget ytarget ztarget]),double([zreference]));
            obj.transformz.forward= findAffineTransformZ(double([xreference  yreference zreference]),double([ztarget]));
        end
        
        function [xo,yo,zo]=transformCoordinatesFwd(obj,x,y,z)
            %transforms reference onto target
%             [xo,yo]=transformPointsForward(obj.transform.forward,x,y);
            x=x/1000;y=y/1000;
            [xo,yo]=transformPointsInverse(obj.transform.inverse,x,y); %inverse of inverse is forward
           
            if nargin>3 %z coordinates present
                z=z/1000;
                 X=transformPointsInverse(obj.transformz.inverse,[x,y,z]);
                 xo=X(:,1);yo=X(:,2);zo=X(:,3);
                 zo=zo*1000;
            else
                zo=[];
            end
             xo=xo*1000;yo=yo*1000;
            
        end
        function [xo,yo,zo]=transformCoordinatesInv(obj,x,y,z)
            %transforms target onto reference
%             [xo,yo]=transformPointsForward(obj.transform.inverse,x,y);
            x=x/1000;y=y/1000;
            [xo,yo]=transformPointsInverse(obj.transform.forward,x,y);
             if nargin>3 %z coordinates present
                z=z/1000;
                 X=transformPointsInverse(obj.transformz.inverse,[x,y,z]);
                 xo=X(:,1);yo=X(:,2);zo=X(:,3);
                 zo=zo*1000;
            else
                zo=[];
            end
            xo=xo*1000;yo=yo*1000;
        end        
        function imout=transformImageFwd(obj,image,cam_pixnm,roi)
            imout=transformImage(obj.transform.forward,image,cam_pixnm,roi);
        end

        function imout=transformImageInv(obj,image,cam_pixnm,roi)
            imout=transformImage(obj.transform.inverse,image,cam_pixnm,roi);
        end  
%         function srim=getSrimage(obj)
%         end
        function ind=getRef(obj,x,y) %returns indices for reference.
            separator=obj.tinfo.separator;
            switch obj.tinfo.targetpos
                case 'top'
                    ind=y>separator(2);
                case 'bottom'
                    ind=y<separator(2);
                case 'left'
                    ind=x>separator(1);
                case 'right'
                    ind=x<separator(1);
            end  
        end
        function makeAffine2d(obj,A)
            tform=affine2d(A);
            obj.transform.forward=tform;
            obj.transform.inverse=invert(tform);
            mirror=struct('midmirror',17664,'targetmirror','no mirror');
            obj.tinfo=struct('targetpos','all','separator',[70656 35328],'mirror',mirror);
        end
    end
end



