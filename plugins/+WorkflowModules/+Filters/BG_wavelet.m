classdef BG_wavelet<interfaces.WorkflowModule
    properties
    end
    methods
        function obj=BG_wavelet(varargin)
            obj@interfaces.WorkflowModule(varargin{:});;
            obj.outputChannels=1; %1: image-background. 2: background image
            obj.outputParameters={'loc_subtractbg'};
        end
        function pard=pardef(obj)
            pard=pardef;
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
            obj.setPar('loc_bg_dt',1); %preview: tiff loader should load only one frame
        end
        function prerun(obj,p)
        end
        function output=run(obj,data,p)
           
            if ~p.loc_subtractbg
                dato2=data;
                dato2.data=(0*data.data);
                output=dato2;
                return
            end
            if ~isempty(data.data)
                img=data.data;
                bg=mywaveletfilter(img,p.loc_wavelet_level);
                
                dato=data;
                dato.data=bg;
                output=dato; 
            else %eof
                output=data;
            end 
        end
       
    end
end


function pard=pardef
pard.loc_subtractbg.object=struct('Style','checkbox','String','Subtract background','Value',1);
pard.loc_subtractbg.position=[1,1];
pard.loc_subtractbg.Width=2;

pard.text1.object=struct('Style','text','String','Wavelet filtering:');
pard.text1.position=[2,1];

pard.text2.object=struct('Style','text','String','Wavelet level');
pard.text2.position=[3,1.3];
pard.loc_wavelet_level.object=struct('Style','edit','String','3');
pard.loc_wavelet_level.position=[3,2.3];
pard.loc_wavelet_level.Width=.7;

end