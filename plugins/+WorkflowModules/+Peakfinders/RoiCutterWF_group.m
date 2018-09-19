classdef RoiCutterWF_group<interfaces.WorkflowModule
    properties
        loc_ROIsize
        preview
        temprois
        tempinfo
        dX
        dT
        dn
    end
    methods
        function obj=RoiCutterWF_group(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=2; 
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
            obj.addSynchronization('loc_ROIsize',obj.guihandles.loc_ROIsize,'String')
            obj.setInputChannels(2,'frame');

        end
        function prerun(obj,p) 
             global tempinfo temprois
            p=obj.getAllParameters;
            obj.loc_ROIsize=p.loc_ROIsize;
            obj.preview=obj.getPar('loc_preview');
            obj.setPar('loc_ROIsize',p.loc_ROIsize);
            obj.dX=p.loc_dTdx(2);
            obj.dT=p.loc_dTdx(1);
            kernelSize=obj.loc_ROIsize;
            obj.dn=ceil((kernelSize-1)/2);
            ninit=100;
            init=zeros(ninit,1);
            tempinfo=struct('inuse',false(size(init)),'x',init,'y',init,'dT',init,'numrois',init,'inframes',{{}});
            temprois=zeros(obj.loc_ROIsize,obj.loc_ROIsize,ninit,'single');
           
        end
        function outputdat=run(obj,data,p)
            persistent tempinfo temprois
            if isempty(tempinfo)
                 ninit=100;
                 init=zeros(ninit,1);
                tempinfo=struct('inuse',false(size(init)),'x',init,'y',init,'dT',init,'numrois',init,'inframes',{{}});
                temprois=zeros(obj.loc_ROIsize,obj.loc_ROIsize,ninit,'single');
            end
            %for challange: fit two times (with/wo this) or include info
            %about frames 'on' and then resdistribute.
            %this is anyways needed e.g. for color assignment
            outputdat=[];
            image=data{1}.data;%get;
            frameh=data{1}.frame;
            if ~isempty(image)     
            maxima=data{2}.data;%get;
            if isempty(maxima.x)
                return;
            end
            

            dX=obj.dX;dT=obj.dT;dn=obj.dn;
            for k=1:length(maxima.x)
                addroi(image,maxima.x(k),maxima.y(k),frameh,dX,dT,dn)
            end 
            
          
            
            if obj.preview 
                outputfig=obj.getPar('loc_outputfig');
                if ~isvalid(outputfig)
                    outputfig=figure(209);
                    obj.setPar('loc_outputfig',outputfig);
                end
                outputfig.Visible='on';
                figure(outputfig)
                hold on
                col=[0.3 0.3 0.3];
                    ax=gca;
                for k=1:length(maxima.x)
                    pos=[maxima.x(k)-obj.dn maxima.y(k)-obj.dn maxima.x(k)+obj.dn maxima.y(k)+obj.dn ];
                    plotrect(ax,pos,col);
                end
                [cutoutimages,maximap]=purgeall();
            else
                [cutoutimages,maximap]=purgerois();
            end 
             info=maximap;
            
            info.frame=maximap.x*0+frameh;
            if ~isempty(cutoutimages)
                outs.info=info;
                outs.img=cutoutimages;
                dato=data{1};%.copy;
                dato.data=outs;%set(outs);
                outputdat=dato; 
            end
            
            else
                outputdat=data{1};
            end
        
        function addroi(image,x,y,frame,dX,dT,dn)
%             global tempinfo temprois
%             tempinfo=obj.tempinfo;
%             temprois=obj.temprois;
            inuse=tempinfo.inuse;
            %later: keep x0 as cut out position (has to be same), but
            %update xsearch
%             dX=obj.dX;
%             inxy=find((tempinfo.x(inuse)-x).^2+(tempinfo.y(inuse)-y).^2<dX^2);
            
            indtemp= find(inuse & tempinfo.x>x-dX & tempinfo.x<x+dX & tempinfo.y>y-dX & tempinfo.y<y+dX,1,'first');
            if ~isempty(indtemp) %already there
%                 finuse=find(inuse);
%                 indtemp=finuse(inxy(1)); %later: choose closest
                tim=cutoutimage(image,tempinfo.x(indtemp),tempinfo.y(indtemp),dn);
                temprois(:,:,indtemp)=temprois(:,:,indtemp)+tim;
                %fill info
                tempinfo.dT(indtemp)=dT; %reset
                tempinfo.numrois(indtemp)=tempinfo.numrois(indtemp)+1;
                tempinfo.inframes{indtemp}(end+1)=frame;
            else %new ROI, this would be standard
                newind=find(tempinfo.inuse==false,1,'first');
                if isempty(newind)
                    newind=length(tempinfo.inuse)+1;
                end
                temprois(:,:,newind)=cutoutimage(image,x,y,dn);
                tempinfo.dT(newind)=dT; %reset
                tempinfo.numrois(newind)=1;
                tempinfo.x(newind)=x;tempinfo.y(newind)=y;  
                tempinfo.inuse(newind)=true;
                tempinfo.inframes{newind}=frame;
            end
        end

        function [cutoutimages,maximap]=purgerois()
%               global tempinfo temprois

            finuse=find(tempinfo.inuse);
            indout=tempinfo.dT(finuse)<0;
            fout=finuse(indout);
            
            cutoutimages=temprois(:,:,fout);
            maximap.x=tempinfo.x(fout);
            maximap.y=tempinfo.y(fout);
            maximap.inframes=tempinfo.inframes(fout);
            tempinfo.inuse(fout)=false;
            tempinfo.dT(finuse)=tempinfo.dT(finuse)-1;% count down dark frames
        end
         function [cutoutimages,maximap]=purgeall()
%              global tempinfo temprois
            finuse=find(tempinfo.inuse);
            fout=finuse;
            cutoutimages=temprois(:,:,fout);
            maximap.x=tempinfo.x(fout);
            maximap.y=tempinfo.y(fout);
            maximap.inframes=tempinfo.inframes(fout);
            
            tempinfo.inuse(fout)=false;
%             tempinfo.dT(finuse)=tempinfo.dT(finuse)-1;% count down dark frames
         end
        end
    end
end

        function out=cutoutimage(image, x,y,dn)
            sim=size(image);
                if  (y<=dn||y>sim(1)-dn||x<=dn||x>sim(2)-dn)
                    y=max(dn+1,min(y,sim(1)-dn));
                    x=max(dn+1,min(x,sim(2)-dn));
                end  
                out=image(y-dn:y+dn,x-dn:x+dn);
        end

function pard=guidef
pard.text.object=struct('Style','text','String','Size ROI (pix)');
pard.text.position=[1,1];


pard.loc_ROIsize.object=struct('Style','edit','String','7');
pard.loc_ROIsize.position=[2,1];
pard.loc_ROIsize.Width=0.7;
pard.loc_ROIsize.TooltipString=sprintf('Size (pixels) of regions around each peak candidate which are used for fitting. \n Depends on fitter. Use larger ROIs for 3D data.');

pard.loc_dTdx.object=struct('Style','edit','String','1 1.5');
pard.loc_dTdx.position=[2,1.7];
pard.loc_dTdx.TooltipString=sprintf('Groupping parameters. dT (frames), dX (pixels)');
pard.loc_dTdx.Width=0.3;

pard.syncParameters={{'loc_ROIsize','loc_ROIsize',{'String'}}};

pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='This plugin cuts out regions of interest of a defined size around the candidate positions and passes these on to the fitter';
end