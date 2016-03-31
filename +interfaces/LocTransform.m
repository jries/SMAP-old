classdef LocTransform<handle
    properties
        transform %forward .inversie
        tinfo
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
        function [xo,yo]=transformCoordinatesFwd(obj,x,y)
            %transforms reference onto target
%             [xo,yo]=transformPointsForward(obj.transform.forward,x,y);
            x=x/1000;y=y/1000;
            [xo,yo]=transformPointsInverse(obj.transform.inverse,x,y); %inverse of inverse is forward
            xo=xo*1000;yo=yo*1000;
            
        end
        function [xo,yo]=transformCoordinatesInv(obj,x,y)
            %transforms target onto reference
%             [xo,yo]=transformPointsForward(obj.transform.inverse,x,y);
            x=x/1000;y=y/1000;
            [xo,yo]=transformPointsInverse(obj.transform.forward,x,y);
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
        
%         function findTransform_out(obj,xreference,yreference,xtarget,ytarget)
%             switch obj.transformType
%                 case {'lwm','polynomial'}
%                     obj.transform = fitgeotrans(double([xtarget ytarget]), [double(xreference) double(yreference)],obj.transformType,obj.transformParameter);
%                 otherwise
%                     obj.transform = fitgeotrans(double([xtarget ytarget]), [double(xreference) double(yreference)],obj.transformType);
%             end
%         end
%         function [xo,yo,ind]= transformCoordinates_out(obj,x,y,part)
%             % part= forward, inverse, left, right, top, bottom
%             if nargin <4
%                 part='inverse';
%             end
%             [ind,direction]=obj.getPart(x,y,part);
%             if direction==-1
%                 trafo=@transformPointsForward;
%             else
%                 trafo=@transformPointsInverse;
%             end
% %              [x,y]=obj.getmirror(x,y,part);
%             [xo,yo]=trafo(obj.transform,x(ind),y(ind));
%         end
%         function imout= transformImages(obj,image,cam_pixnm,part)
%             %size(imout)=size(image)
%             %part = reference, target
% %             image=image';
% %             disp('not implemented yet')
%             %not use imref object: change transformation using cam_pixnm
%             %imwarp with XData and YData
%             
%             sizeim=size(image);
%             R = imref2d(sizeim,cam_pixnm,cam_pixnm);
%             imout = imwarp(image,R,obj.transform,'OutputView',R);
% %             imout=imwarp(image,obj.transform,'bicubic','XData',[0 sizeim(1)],'YData',[0 sizeim(2)]);
%             
%         end
%         function range=getRange(obj,part,roi)
%             switch part
%                 case {'forward','inverse','all'} %all
%                     range.x=[roi(1) roi(1)+roi(3)];
%                     range.y=[roi(2) roi(2)+roi(4)];
%                 case 'reference' %left
%                     switch obj.refPart
%                         case 'all'        
%                             range.x=[roi(1) roi(1)+roi(3)];
%                             range.y=[roi(2) roi(2)+roi(4)];
%                         case 'left-right'
%                             range.x=[roi(1) obj.midpoint];
%                             range.y=[roi(2) roi(2)+roi(4)];
%                         case 'top-bottom';
%                             range.x=[roi(1) roi(1)+roi(3)];
%                             range.y=[roi(2) obj.midpoint];
%                     end
%                 case 'target' %right
%                     switch obj.refPart
%                         case 'all'        
%                             range.x=[roi(1) roi(1)+roi(3)];
%                             range.y=[roi(2) roi(2)+roi(4)];
%                         case 'left-right'
%                             range.x=[obj.midpoint obj.midpoint+(obj.midpoint-roi(1))];
%                             range.y=[roi(2) roi(2)+roi(4)];
%                         case 'top-bottom';
%                             range.x=[roi(1) roi(1)+roi(3)];
%                             range.y=[obj.midpoint obj.midpoint+(obj.midpoint-roi(2))];
%                     end                  
%             end
%         end
%         function [locsout,ind]=getPartLocs(obj,locs,x,y,part)
%             
%             ind=obj.getPart(x,y,part);
%             fn=fieldnames(locs);
%             for k=1:length(fn)
%                 locsout.(fn{k})=locs.(fn{k})(ind);
%             end
%             [locsout.x,locsout.y]=obj.getmirror(locsout.x,locsout.y,part);
%         end
%         function [ind,direction]=getPart(obj,x,y,part)
%             %new: refpart: up-down or right-left or all
%             % part can be forward, backward, reference, target
%             switch part
%                 case 'forward' %all
%                     ind=true(size(x));
%                     direction=1;
%                 case 'inverse'
%                     ind=true(size(x));
%                     direction=-1;
%                 case 'reference' %left
% %                     ind=x<=obj.midpoint;
%                     direction=1;
%                     switch obj.refPart
%                         case 'all'        
%                             ind=true(size(x));
%                         case 'left-right'
%                             ind=x<=obj.midpoint;
%                         case 'top-bottom';
%                             ind=y<=obj.midpoint;
%                     end
%                 case 'target' %right
%                     direction =-1;
%                     switch obj.refPart
%                         case 'all'        
%                             ind=true(size(x));
%                         case 'left-right'
%                             ind=x>obj.midpoint;
%                         case 'top-bottom';
%                             ind=y>obj.midpoint;
%                     end                   
%             end
%         end
%         function [srim,range]=getSRimage(obj,xall,yall,indin,part,pixelsize,roi)
%             indpart=obj.getPart(xall,yall,part);
%             
%             x=xall(indpart&indin);
%             y=yall(indpart&indin);
%             [x,y]=obj.getmirror(x,y,part);
%             range=obj.getRange(part,roi);
%             srim=myhist2(x,y,pixelsize,pixelsize,range.x,range.y);
%          
% 
%         end
%         function [x,y]=getmirror(obj,x,y,part)
%         if strcmp(part,'target')
%             switch obj.mirror.direction
%                 case {'x' }
%                     x=2*obj.mirror.midpoint(1)-x;
%                     
%                 case {'y'} 
%                     y=2*obj.mirror.midpoint(2)-y;
%                     
%                 case 'xy'
%                     x=2*obj.mirror.midpoint(1)-x;
%                     y=2*obj.mirror.midpoint(2)-y;
%                 otherwise
%                     disp('none')
%                   
%             end
% 
%         end
            
%         end
        

    end
end


function imout=transformImage(transformf,image,cam_pixnm,roi)
    if nargin<4
        roi=[0 0 size(image)];
    end
    sizeim=size(image);
    extxnm=[roi(1) roi(1)+roi(3)]*cam_pixnm/1000;
    extynm=[roi(2) roi(2)+roi(4)]*cam_pixnm/1000;
%             R = imref2d(sizeim,cam_pixnm,cam_pixnm);
    R = imref2d(sizeim,extxnm,extynm);
    imout = imwarp(image,R,transformf,'OutputView',R);
end
