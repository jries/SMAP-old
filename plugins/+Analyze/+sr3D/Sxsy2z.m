classdef Sxsy2z<interfaces.DialogProcessor
    properties
        cal3D
        outsx2sy2
    end
    methods
        function obj=Sxsy2z(varargin)        
           obj@interfaces.DialogProcessor(varargin{:}) ;
             obj.inputParameters={'cam_pixelsize_nm'};
        end
        function initGui(obj)
            set(obj.guihandles.savebutton,'Callback',{@loadcalfile,obj})
        end
        function out=run(obj,p)
            out=[];
            locs=obj.locData.getloc({'PSFxnm','PSFynm'});
            if isempty(locs.PSFynm)
                error('PSFynm required')
            end
            sx2sy2=(locs.PSFxnm/p.cam_pixelsize_nm).^2-(locs.PSFynm/p.cam_pixelsize_nm).^2;
            z=feval(obj.outsx2sy2,sx2sy2);
%             refractiveIndexMismatch=1.25
            obj.locData.loc.znm=-z*1000*p.refractiveIndexMismatch;
            obj.locData.regroup;
            obj.setPar('locFields',fieldnames(obj.locData.loc))
            
        end
        function pard=pardef(obj)
            pard=pardef;
        end
    end
end


function loadcalfile(a,b,obj)
fn=obj.guihandles.calfile.String;
[f,p]=uigetfile(fn);
if f
load([p f])
obj.cal3D=cal3D;
obj.outsx2sy2=outsx2sy2;
obj.guihandles.calfile.String=[p f];
end
end

function pard=pardef
% pard.d3_color.object=struct('Style','checkbox','String','render in color');
% pard.d3_color.position=[2,1];

pard.zposobjective_check.object=struct('String','z-position (um) of objective','Style','checkbox');
pard.zposobjective_check.position=[1,1];
pard.zposobjective_check.Width=2;


pard.zposobjective.object=struct('Style','edit','String','30'); 
pard.zposobjective.position=[1,3];

pard.pol_check.object=struct('String','give polynome here:','Style','checkbox');
pard.pol_check.position=[2,1];
pard.pol_check.Width=2;


pard.fitpol.object=struct('Style','edit','String','1 0 0'); 
pard.fitpol.position=[2,3];


pard.t1.object=struct('Style','edit','String','refractive Index Mismatch factor (<=1)'); 
pard.t1.position=[3,1];
pard.refractiveIndexMismatch.object=struct('Style','edit','String','1'); 
pard.refractiveIndexMismatch.position=[3,2];

pard.calfile.object=struct('Style','edit','String','settings/cal_3DAcal.mat');
pard.calfile.position=[4,1];
pard.calfile.Width=3;


pard.savebutton.object=struct('String','load','Style','pushbutton');
pard.savebutton.position=[4,4];



end
