classdef RegisterImages<interfaces.DialogProcessor
    properties
        transformation
    end
    methods
        function obj=RegisterImages(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'sr_pixrec'};
        end
        
        function out=run(obj,p)
            out=[];
            pixrec=p.sr_pixrec;
            targetlayer=p.targetselect.Value;
            reflayer=p.refselect.Value;
            imref=sum(obj.locData.layer(reflayer).images.finalImages.image,3);
            imtarget=sum(obj.locData.layer(targetlayer).images.finalImages.image,3);
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
        function pard=pardef(obj)
            pard=pardef(obj);
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




function pard=pardef(obj)
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
% pard.targetpart.object=struct('Style','popupmenu','String','all|left|right|top|bottom');
% pard.targetpart.position=[3,1];
% pard.refpart.object=struct('Style','popupmenu','String','all|left|right|top|bottom');
% pard.refpart.position=[3,2];

% pard.targetmirror.object=struct('Style','popupmenu','String','no mirror|left-right|up-down');
% pard.targetmirror.position=[4,1];
% pard.refmirror.object=struct('Style','popupmenu','String','no mirror|left-right|up-down');
% pard.refmirror.position=[4,2];

pard.Tfile.object=struct('Style','edit','String','settings/temp/temp_T.mat');
pard.Tfile.position=[8,1];
pard.Tfile.Width=3;

pard.savebutton.object=struct('Style','pushbutton','String','save T','Callback',@obj.savebutton);
pard.savebutton.position=[8,4];
% 
% pard.registeronlocs.object=struct('Style','checkbox','String','register locs');
% pard.registeronlocs.position=[1,3.5];

% pard.transform.object=struct('Style','popupmenu','String','translate|affine|LWM');
% pard.transform.position=[2,3.5];
end