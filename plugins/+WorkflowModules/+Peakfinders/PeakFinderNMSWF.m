classdef PeakFinderNMSWF<interfaces.WorkflowModule
    properties
        NMScutoff
    end
    methods
        function obj=PeakFinderNMSWF(varargin)
            obj@interfaces.WorkflowModule(varargin{:});
            obj.inputParameters={'loc_loc_filter_sigma','EMexcessNoise'};
%             obj.parameters.EMexcessNoise=1;
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
            obj.guihandles.cutoffmode.Callback={@cutoffmode_callback,obj};

             
%            obj.guihandles.camparbutton.Callback={@camparbutton_callback,obj};
        end
        function prerun(obj,p)
            p=obj.getAllParameters;
            filtersize=p.loc_loc_filter_sigma;
            excess=p.EMexcessNoise;
%             if isempty(escess)
%                 excess=1;
%             end
            obj.NMSkernel=ceil((p.NMS_kernel_size-1)/2);
            state=p.cutoffmode;
            val=(p.NMS_cutoff);
            if ~state
                val=sqrt(2)*erfinv(1-2*val); %in paper it is 2p-1 
                PSFx0=1;
                val=val*sqrt(PSFx0^2/(PSFx0^2+filtersize^2)*2*excess);
            end
            obj.NMScutoff=val;
            %  co=co*sqrt(p.PSFx0^2/(p.PSFx0^2+p.filtersize^2)*2);


        end
        function dato=run(obj,data,p)
%             output=[];
%             data.frame
            image=data.data;%get;
            maxima=NMS2DBlockCcall(image,obj.NMSkernel); %find maxima

            maxind= (maxima(:,3)>obj.NMScutoff);
            maxout.y=maxima(maxind,1);
            maxout.x=maxima(maxind,2);
            dato=data;%.copy;
            
            dato.data=(maxout);
        end
    end
end

function val = prob2photon(p,PSFx0,filtersize,excess)
val=sqrt(2)*erfinv(1-2*p); %in paper it is 2p-1 
% PSFx0=1;
 val=val*sqrt(PSFx0^2/(PSFx0^2+filtersize^2)*2*excess);
end

function p=photon2prob(val,PSFx0,filtersize,excess)
val=val/sqrt(PSFx0^2/(PSFx0^2+filtersize^2)*2*excess);
p=(1-erf(val/sqrt(2)))/2;
end

function cutoffmode_callback(a,b,obj)
PSFx0=1;
p=obj.getAllParameters;
filtersize=p.loc_loc_filter_sigma;
excess=p.EMexcessNoise;
state=p.cutoffmode;
 val=p.NMS_cutoff;
if state
    obj.guihandles.cutoffmode.String='Absolute';
    
    p=prob2photon(val,PSFx0,filtersize,excess);
   
else
    obj.guihandles.cutoffmode.String='probability (p<1)';
        p=photon2prob(val,PSFx0,filtersize,excess);

end

p=max(1E-7,p);
p=min(1E10,p);
 obj.guihandles.NMS_cutoff.String=num2str(p);
end

function pard=guidef
pard.cutoffmode.object=struct('Style','togglebutton','String','probability (p<1)');
pard.cutoffmode.position=[1,1];
pard.cutoffmode.Width=1.3;
pard.NMS_cutoff.object=struct('Style','edit','String','0.01');
pard.NMS_cutoff.position=[1,2.3];
pard.NMS_cutoff.Width=.7;
pard.text.object=struct('Style','text','String','NMS kernel size (pix)');
pard.text.position=[2,1];
pard.text.Width=1.3;
pard.NMS_kernel_size.object=struct('Style','edit','String','7');
pard.NMS_kernel_size.position=[2,2.3];
pard.NMS_kernel_size.Width=.7;
end