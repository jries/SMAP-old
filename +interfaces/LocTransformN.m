classdef LocTransformN<handle
    properties
        transform2Reference
        transform2Target
        tinfo
        transformZ2Reference
        transformZ2Target
        
    end
    
    methods
        function obj=LocTransformN(varargin)
            obj.transform2Reference{1}=affine2d(eye(3)); %initialize channel 1 to reference
            obj.transform2Target{1}=affine2d(eye(3));
            obj.transformZ2Reference{1}=affine3d(eye(4)); %initialize channel 1 to reference
            obj.transformZ2Target{1}=affine3d(eye(4));
            tinfo{1}.xrange=[-inf inf];
            tinfo{1}.yrange=[-inf inf];
        end
        function setTransform(obj,varargin)
            properties={'xrange','yrange','type','parameter','mirror'};
            for k=1:2:length(varargin)
                if any(strcmp(properties,varargin{k}))
                    obj.tinfo.(varargin{k})=varargin{k+1};
                else 
                    disp([varargin{k} ' is not a proper parameter']);
                end
            end
        end
        function findTransform(obj,channel,coordreference,coordtarget,type,parameter)
            %XXX make sure reference and target are not exchanged 
            %matrix becomes singular for too large coordinates. Use um
            %instead of nm
            
            if nargin<5
                th=obj.tinfo{channel};
                if isfield(th,'type') && ~isempty(th.type)
                    type=th.type;
                    if isfield(th,'parameter')
                        parameter=th.parameter;
                    end
                else
                    parameter=[];
                    for k=1:length(obj.tinfo) %look for other channels if any type is specified, use the same
                         th=obj.tinfo{k};
                        if isfield(th,'type') && ~isempty(th.type)
                            type=th.type;
                            if isfield(th,'parameter')
                                parameter=th.parameter;
                            end
                            break
                        end
                        
                    end
                    if isempty(parameter)
                        warning('no transformation type specified');
                        return
                    end
                end
            end
            coordreference=double(coordreference/1000);
            coordtarget=double(coordtarget/1000);
            switch type
                case {'lwm','polynomial'}
                    obj.transform2Target{channel} = fitgeotrans(coordtarget(:,1:2),coordreference(:,1:2),type,parameter);
                    obj.transform2Reference{channel} = fitgeotrans( coordreference(:,1:2),coordtarget(:,1:2),type,parameter);
                otherwise             
                    obj.transform2Target{channel}= fitgeotrans(coordtarget(:,1:2),coordreference(:,1:2),type);
                    obj.transform2Reference{channel}= fitgeotrans(coordreference(:,1:2),coordtarget(:,1:2),type);
            end
            if size(coordreference,2)>2 %3D data set: also do z-transform
                coordtargettransformed=horzcat(obj.transformToReference(channel,coordtarget(:,1:2)),coordtarget(:,3));
                    %only affine3d possible      
                obj.transformZ2Target{channel}= findAffineTransformZ(coordtargettransformed,coordreference(:,3));
                obj.transformZ2Reference{channel}= findAffineTransformZ(coordreference,coordtarget(:,3));
            end
        end
        
        function co=transformToReference(obj,channel,ci)
            ci=ci/1000;
            co=transformPointsInverse(obj.transform2Reference{channel},ci(:,1:2)); %inverse of inverse is forward          
            if size(ci,2)>2 %z coordinates present               
                 X=transformPointsInverse(obj.transformZ2Reference{channel},horzcat(co,ci(:,3)));
                 co(:,3)=X(:,3);
            end
             co=co*1000; %back to nm   
        end
        function co=transformToTarget(obj,channel,ci)
            ci=ci/1000;
            co=transformPointsInverse(obj.transform2Target{channel},ci(:,1:2)); %inverse of inverse is forward          
            if size(ci,2)>2 %z coordinates present               
                 X=transformPointsInverse(obj.transformZ2Target{channel},horzcat(co,ci(:,3)));
                 co(:,3)=X(:,3);
            end
             co=co*1000; %back to nm   
        end  

        function imout=transformImageToTarget(obj,channel,image,cam_pixnm,roi)
            imout=transformImage(obj.transform2Target{channel},image,cam_pixnm,roi);
        end

        function imout=transformImageToReference(obj,image,cam_pixnm,roi)
            imout=transformImage(obj.transform2Reference{channel},image,cam_pixnm,roi);
        end  
        
        function ind=getPart(obj,channel,coordinates)
            th=obj.tinfo{channel};
            if ~isfield(th,'xrange')||isempty(th.xrange)||~isfield(th,'yrange')||isempty(th.yrange)
                disp('no range specified, getPart in LocTransformN returns all coordinates');
                ind=true(size(coordinates,1),1);
            else
                ind=coordinates(:,1)>th.xrange(1) & coordinates(:,1)<=th.xrange(2) & coordinates(:,2)>th.yrange(1)& coordinates(:,2)<th.yrange(2); 
            end
        end


        function makeAffine2d(obj,channel,A)
            tform=affine2d(A);
            obj.transform2Target=tform;
            obj.transform2Reference=invert(tform);
        end
    end
end



