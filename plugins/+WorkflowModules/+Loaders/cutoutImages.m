classdef cutoutImages<interfaces.WorkflowModule
    properties
        settings3D
    end
    methods
        function obj=cutoutImages(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function prerun(obj,p)
            if exist(p.cal_3Dfile,'file')
                [~,~,ext]=fileparts(p.cal_3Dfile);
                switch ext
                    case '.mat'
                        l=load(p.cal_3Dfile);
                        obj.settings3D=l.parameters.settings_3D;
                    case '.txt'
                        obj.settings3D=readstruct(p.cal_3Dfile);
                    otherwise
                        disp('cannot recognize format')
                end
            else
                obj.settings3D=[];
            end
        end
        function datao=run(obj,data,p)
           img=data.data;
           datao=data;
           if ~isempty(obj.settings3D) && ~isempty(img)
               datao.data=cutoutchannels(img,obj.settings3D);
           end
        end
        function loadcalfile(obj,a,b)
            file=obj.getSingleGuiParameter('cal_3Dfile');
            if ~contains(file,'txt') && ~contains(file,'mat')
                file='*.*';
            end
            [file path] =uigetfile(file,'select 3D settings file for cutout');
            if file
                obj.setGuiParameters(struct('cal_3Dfile',[path file]));
                obj.setPar('cal_3Dfile',[path file]);
            end
        end
    end
end



function pard=guidef(obj)
pard.loadcal.object=struct('Style','pushbutton','String','Load 3D settings','Callback',@obj.loadcalfile);
pard.loadcal.position=[2,1];
pard.loadcal.Width=1;
pard.cal_3Dfile.object=struct('Style','edit','String','');
pard.cal_3Dfile.position=[2,2];
pard.cal_3Dfile.Width=3;
pard.cal_3Dfile.TooltipString=sprintf('settings_3D.txt file specifying ROIs on the chip and mirroring, or 3Dcal.mat file containing this.');

pard.syncParameters={{'cal_3Dfile','cal_3Dfile',{'String'}}};

pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='Allows editing metadata, or loading from a file.';
end