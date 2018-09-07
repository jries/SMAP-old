classdef IntLoc2posN<interfaces.WorkflowModule
    properties
        filestruc;
        locs
        roi
    end
    methods
        function obj=IntLoc2posN(varargin)
            obj@interfaces.WorkflowModule(varargin{:});
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
        end
        function prerun(obj,data,p)
            global intLoc2pos_ind2 intLoc2pos_locframes;
            intLoc2pos_ind2=1;
            
            p=obj.getAllParameters;
          
            transform=loadtransformation(obj,p.Tfile);
            obj.locs=obj.locData.getloc({'frame','xnm','ynm','znm','PSFxnm','bg'});
            cpix=obj.locData.files.file(1).info.cam_pixelsize_um*1000;
            x=double(obj.locs.xnm)/cpix(1);  %in pixels
            y=double(obj.locs.ynm)/cpix(end);
            
%             obj.locs.xA=zeros(size(x));obj.locs.xB=zeros(size(x));obj.locs.yA=zeros(size(x));obj.locs.yB=zeros(size(x));
            

            switch p.Tmode.selection
                case 'target'
                    pos=transform.transformToTarget(2,horzcat(x,y));
                    obj.locs.xA=pos(:,1);obj.locs.yA=pos(:,2);
                case 'reference'
                     obj.locs.xA=x;obj.locs.yA=y;
                case 'both'
                    disp('not implemented')

            end

            
            if ~isempty(obj.locs.znm)
          
                    d=0.42;
                    g=-0.2;
                    sx0=1.1;
                    obj.locs.PSFxpix=sigmafromz_simple(obj.locs.znm/1000,[d -g sx0]);
                    obj.locs.PSFypix=sigmafromz_simple(obj.locs.znm/1000,[d g sx0]);
     
            else
                obj.locs.PSFxpix=double(obj.locs.PSFxnm/transform.cam_pixnm{1}(1));
                obj.locs.PSFypix=obj.locs.PSFxpix;
            end
            intLoc2pos_locframes=obj.locs.frame;
            obj.roi=obj.getPar('loc_fileinfo').roi;

        end
        function datout=run(obj,data,p)
            global intLoc2pos_ind2 intLoc2pos_locframes

                if ~data.eof
            lf=length(obj.locs.xA);
            frame=data.frame;
                %find indices for same frame
                ind1=intLoc2pos_ind2;
                while ind1>0&&intLoc2pos_locframes(ind1)<frame && ind1<lf
                    ind1=ind1+1;
                end
                ind1=min(ind1,lf);
                intLoc2pos_ind2=ind1;
                if ~(intLoc2pos_locframes(intLoc2pos_ind2)==frame) %no localizatiaon in frame
                      datout=data;%.copy;
                     datout.data=struct('x',[]);%.set(maxout);
                    return
                end
                while intLoc2pos_ind2<=lf&&intLoc2pos_locframes(intLoc2pos_ind2)==frame
                    intLoc2pos_ind2=intLoc2pos_ind2+1;
                end
                intLoc2pos_ind2=intLoc2pos_ind2-1;
                if (intLoc2pos_ind2-ind1+1)~=(sum(obj.locs.frame==frame))
                    disp((intLoc2pos_ind2-ind1+1)-(sum(obj.locs.frame==frame)))
                    disp(frame)
                end
                xrel=(obj.locs.xA(ind1:intLoc2pos_ind2)-obj.roi(1));
               yrel=(obj.locs.yA(ind1:intLoc2pos_ind2)-obj.roi(2));
             
                
               maxout.x=round(xrel);
               maxout.y=round(yrel);
               maxout.frame=frame+0*maxout.y;
               maxout.dx=xrel-round(xrel);
               maxout.dy=yrel-round(yrel);
               maxout.PSFxpix=obj.locs.PSFxpix(ind1:intLoc2pos_ind2);
               maxout.PSFypix=(obj.locs.PSFypix(ind1:intLoc2pos_ind2));
               datout=data;%.copy;
               datout.data=maxout;%.set(maxout);
%                obj.output(datout); 
                else
             
                   datout=data;
                    datout.data=struct('x',[]);
                    datout.eof=true;
                end
                

        end
    end
end

% function [loc,locr]=nm2pixLoc(x,y,pixelsize,roi)
% loc.x=(x/pixelsize(1))-roi(1);
% loc.y=(y/pixelsize(2))-roi(2);
% locr.x=round(loc.x);
% locr.y=round(loc.y);
% end


function pard=guidef
pard.Tfile.object=struct('Style','edit','String','Tfile');
pard.Tfile.position=[1,1];
pard.Tfile.Width=1.3;
pard.Tmode.object=struct('Style','popupmenu','String',{{'target','reference','both'}});
pard.Tmode.position=[1,1.3];
pard.Tmode.Width=.7;
% pard.pixcam.object=struct('Style','edit','String','138');
% pard.pixcam.position=[2,1];
% pard.PSFxnm.object=struct('Style','edit','String','138');
% pard.PSFxnm.position=[2,1];
pard.plugininfo.type='WorkflowModule'; 
end

function PSFx=sigmafromz_simple(z,p)%[d g sx0]);
    PSFx=p(3).*sqrt(1+(z-p(2)).^2./p(1).^2);
end


function s=sigmafromz(par,z,B0)
par=real(par);
% parx= [d sx0 Ax Bx g mp]
s0=par(2);d=par(1);A=par(3);B=par(4)*B0;g=par(5);mp=par(6);

% s=s0*sqrt(1+(z-g+mp).^2/d^2);
s=s0*sqrt(1+(z-g+mp).^2/d^2+A*(z-g+mp).^3/d^3+B*(z-g+mp).^4/d^4);
s=real(s);
end