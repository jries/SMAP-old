classdef RegisterImages<interfaces.DialogProcessor
    properties
        transformation
    end
    methods
        function obj=RegisterImages(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'sr_pixrec'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)
            out=[];
            pixrec=p.sr_pixrec;
            targetlayer=p.targetselect.Value;
            reflayer=p.refselect.Value;
            ll=length(obj.locData.layer);
            if targetlayer>ll||reflayer>ll
                out.error=('selected layer does not exist or has not been calculated');
                return;
            end
            imref=sum(obj.locData.layer(reflayer).images.finalImages.image,3);
            imtarget=sum(obj.locData.layer(targetlayer).images.finalImages.image,3);
            if sum(imref(:))==0||sum(imtarget(:))==0
                out.error=('reconstructed images are empty.');
                return;
            end
            initaxis(p.resultstabgroup,'correlation');
             [dy,dx]=getShiftCorr(imref,imtarget,1);
             dx=dx*pixrec;
             dy=dy*pixrec;
            dxt=obj.getPar('shiftxy_min','layer',targetlayer);
            dyt=obj.getPar('shiftxy_max','layer',targetlayer);
%             dxr=obj.getPar('shiftxy_min','layer',reflayer);
%             dyr=obj.getPar('shiftxy_max','layer',reflayer);
            obj.setPar('shiftxy_min','layer',targetlayer,dx+dxt(1))
            obj.setPar('shiftxy_max','layer',targetlayer,dy+dyt)
            notify(obj.P,'sr_render');
            
            %make locTransform
            A=[1,0,0;0,1,0;(dx+dxt(1))/1000,(dy+dyt)/1000,1];
            trafo=interfaces.LocTransform;
            trafo.makeAffine2d(A);
            obj.transformation=trafo;
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function initGui(obj)
            obj.addSynchronization('transformationfile',obj.guihandles.Tfile,'String');
        end
        function savebutton(obj,a,b)
            fn=obj.guihandles.Tfile.String;
            [f,path]=uiputfile(fn,'Select transformation file _T.mat');
            if f
                obj.guihandles.Tfile.String=[path f];
                transformation=obj.transformation;
                save([path,f],'transformation');
            end      
        end 
    end
end




function pard=guidef(obj)
pard.texta.object=struct('String','target','Style','text');
pard.texta.position=[1,1];
pard.textb.object=struct('String','reference','Style','text');
pard.textb.position=[1,2];
pard.targetselect.object=struct('Style','popupmenu','String','layer1|layer2|layer3|layer4|layer5');
pard.targetselect.position=[2,1];
pard.targetselect.load=false;
pard.refselect.object=struct('Style','popupmenu','String','layer1|layer2|layer3|layer4|layer5');
pard.refselect.position=[2,2];
pard.refselect.load=false;

pard.Tfile.object=struct('Style','edit','String','settings/temp/temp_T.mat');
pard.Tfile.position=[8,1];
pard.Tfile.Width=3;

pard.savebutton.object=struct('Style','pushbutton','String','save T','Callback',@obj.savebutton);
pard.savebutton.position=[8,4];

pard.plugininfo.description=sprintf(['Register Images calculates shift between rendered images in two layers and writes this shift into the channel tab of the target layer.'...
    '\n Adjust pixel size and FoV so that the shift can be calculated from the reconstructed image.'...
    '\n Transformation can also be saved for later use with Apply Transform']);
pard.plugininfo.name='Register Images';
end