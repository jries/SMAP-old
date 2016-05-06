classdef ImageFilter<interfaces.WorkflowModule
    properties
        filterkernel
        preview
    end
    methods
        function obj=ImageFilter(varargin)
            obj@interfaces.WorkflowModule(varargin{:});
            obj.outputParameters={'loc_loc_filter_sigma'};
        end
        function pard=guidef(obj)
            pard=guidef;
        end
%         function initGui(obj)
%             initGui@interfaces.WorkflowModule(obj);
% %            obj.guihandles.loadmetadata.Callback={@loadmetadata_callback,obj};
% %            obj.guihandles.camparbutton.Callback={@camparbutton_callback,obj};
%         end
        function prerun(obj,p)
            p=obj.getAllParameters;
            fs=p.loc_loc_filter_sigma;
            if fs>0
            obj.filterkernel=fspecial('gaussian', max(5,ceil(3.5/2*fs)*2+1), fs);
            else
                obj.filterkernel=1;
            end
            obj.preview=obj.getPar('loc_preview');
        end
        function dato=run(obj,data,p)
            dato=data;%.copy;
            imf=filter2(obj.filterkernel,data.data);
            if obj.preview&&obj.getPar('loc_previewmode').Value==3&&~isempty(imf)
              
                figure(obj.getPar('loc_outputfig'))
                hold off
                imagesc(imf);
                colorbar;
                axis equal
            end
            dato.data=(imf);
        end
    end
end


function pard=guidef
pard.text.object=struct('Style','text','String','Filter: Sigma (pix): ');
pard.text.position=[1,1];
pard.text.Width=1.3;
pard.loc_loc_filter_sigma.object=struct('Style','edit','String','1.2');
pard.loc_loc_filter_sigma.position=[1,2.3];
pard.loc_loc_filter_sigma.Width=.7;
pard.loc_loc_filter_sigma.TooltipString=sprintf('Sigma (in camera pixels) for a Gaussian filter which is applied after \n background correction and before peak finding.');
pard.plugininfo.type='WorkflowModule';
end