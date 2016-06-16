classdef RoiCutterWF<interfaces.WorkflowModule
    properties
        loc_ROIsize
        preview

    end
    methods
        function obj=RoiCutterWF(varargin)
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
            
            p=obj.getAllParameters;
            obj.loc_ROIsize=p.loc_ROIsize;
            obj.preview=obj.getPar('loc_preview');
           
        end
        function outputdat=run(obj,data,p)
            outputdat=[];
            image=data{1}.data;%get;
            if ~isempty(image)     
            maxima=data{2}.data;%get;
            if isempty(maxima.x)
                return;
            end
%             data{1}.frame
            
            
            kernelSize=obj.loc_ROIsize;
            dn=ceil((kernelSize-1)/2);
            sim=size(image);
            cutoutimages=zeros(kernelSize,kernelSize,length(maxima.x),'single');
          
            
%             info.x=zeros(length(maxima.x),1,'single');
%             info.y=zeros(length(maxima.x),1,'single');
%             info.frame=zeros(length(maxima.x),1,'single');
            ind=0;
            goodind=~(maxima.y<=dn|maxima.y>sim(1)-dn|maxima.x<=dn|maxima.x>sim(2)-dn);
%             info(sum(goodind))=struct('x',[],'y',[],'frame',[]);
%             info(length(maxima.x))=struct('x',[],'y',[],'frame',[]);
            for k=1:length(maxima.x)
                 ind=ind+1;
                if goodind(k)
%                 if maxima.y(k)>dn&&maxima.y(k)<=sim(1)-dn&&maxima.x(k)>dn&&maxima.x(k)<=sim(2)-dn;
                   
%                     ind=k; %return empty frame if outside
                    cutoutimages(:,:,ind)=image(maxima.y(k)-dn:maxima.y(k)+dn,maxima.x(k)-dn:maxima.x(k)+dn); %coordinates exchanged.
                else
                    cutoutimages(:,:,ind)=zeros(kernelSize,kernelSize,1,'single');
                end
%                     info(ind).x=maxima.x(k);
%                     info(ind).y=maxima.y(k);
%                     info(ind).frame=data{1}.frame;     
            end 
            info.x=maxima.x;
            info.y=maxima.y;
            frameh=data{1}.frame;
            info.frame=maxima.x*0+frameh;
           

            outs.info=info;
            outs.img=cutoutimages(:,:,1:ind);
            dato=data{1};%.copy;
            dato.data=outs;%set(outs);
%             obj.output(dato)
            outputdat=dato;

            
            
            
            if obj.preview && ~isempty(maxima.x)
%                 figure(obj.globpar.parameters.outputfig)
%                 ax=gca;
                ax=findobj(obj.getPar('loc_outputfig').Children,'Type','Axes');
                ax.NextPlot='add';
%                 hold on
                for k=1:length(maxima.x)
                    pos=[maxima.x(k)-dn maxima.y(k)-dn maxima.x(k)+dn maxima.y(k)+dn ];
                    col=[0.3 0.3 0.3];
                    plotrect(ax,pos,col);
                end

%             drawnow
            
            end 
            else
%                 obj.output(data{1});
                outputdat=data{1};
            end
%             outputdat=[];
        end
        

    end
end



function pard=guidef
pard.text.object=struct('Style','text','String','Size ROI (pix)');
pard.text.position=[1,1];


pard.loc_ROIsize.object=struct('Style','edit','String','7');
pard.loc_ROIsize.position=[2,1];
pard.loc_ROIsize.TooltipString=sprintf('Size (pixels) of regions around each peak candidate which are used for fitting. \n Depends on fitter. Use larger ROIs for 3D data.');
pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='This plugin cuts out regions of interest of a defined size around the candidate positions and passes these on to the fitter';
end