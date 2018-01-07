classdef applyexistingdriftcorrection<interfaces.DialogProcessor
    methods
        function obj=applyexistingdriftcorrection(varargin)        
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.history=true;
                obj.showresults=false;
                obj.guiselector.show=true;
        end
        
        function out=run(obj,p)
            out=[];
            obj.setPar('undoModule','applydriftcorrection');
            notify(obj.P,'backup4undo');
            try
                saveloc=load(p.file);
                driftinfo=saveloc.saveloc.file(1).driftinfo; 
            catch err
                disp('could not read file or no drift info present');
                return
            end

            if p.correctxy&&isfield(driftinfo,'xy')
                dr.xy=driftinfo.xy;
                dr.xy.x=driftinfo.xy.dxt;
                dr.xy.y=driftinfo.xy.dyt;
            end
            if p.correctz&&isfield(driftinfo,'z')
                dr.z=driftinfo.z;
                dr.z.z=driftinfo.z.zt;
            end
            locsnew=applydriftcorrection(dr,obj.locData.loc);
            obj.locData.loc.xnm=locsnew.xnm;
            obj.locData.loc.ynm=locsnew.ynm;

            obj.locData.files.file(1).driftinfo=driftinfo;
            fn=obj.locData.files.file(1).name;
            if strfind(fn,'_sml')
                fnn=strrep(fn,'_sml','_adriftc_sml');
            else
                fnn=strrep(fn,'fitpos','adriftc_sml');
            end
            if p.save_dc
                obj.locData.savelocs(fnn); 
            end
            obj.locData.regroup;
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function loadfile(a,b,obj)
oldf=obj.getSingleGuiParameter('file');
if ~exist(oldf,'file')
    oldf=obj.getPar('lastSMLFile');
end
[f,path]=uigetfile(oldf);
if f
    obj.setGuiParameters(struct('file',[path f]));
end
end


function pard=guidef(obj)

pard.correctxy.object=struct('String','Correct xy-drift','Style','checkbox','Value',1);
pard.correctxy.position=[1,1];

pard.correctz.object=struct('String','Correct z-drift','Style','checkbox','Value',0);
pard.correctz.position=[2,1];


pard.file.object=struct('String','','Style','edit');
pard.file.position=[3,1];
pard.file.Width=3.5;

pard.loadfile.object=struct('String','load','Style','pushbutton','Callback',{{@loadfile,obj}});
pard.loadfile.position=[3,4.5];
pard.loadfile.Width=0.5;

pard.save_dc.object=struct('String','Save driftcorrected SML','Style','checkbox','Value',0);
pard.save_dc.position=[4,1];
pard.save_dc.Width=2;
pard.save_dc.Optional=true;

pard.plugininfo.name='apply drift correction from driftcorrected file';
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description={'Drift correction based on cross-correlation.','Algorithm: the data set is divided into [timepoints] blocks, for which superresolution images are calculated. The displacement between all images is calcualted with a FFT-based cross-correlation algorithm. The position of the maxima of the cross-correlation curve are fitted with sub-pixel accuracy with a free elliptical Gaussian.',...
    'A robust estimator is used to calculate the drift vs frame from all pairwise displacements.','All localiaztions visible in the superresolution image are used to infer the drift. Use [Render]...[Layer] to control this.',...
    'If two files are loaded, their drift is calculated together and they are saved as one file with their filenumbers copied to the channel field.',' ','(c) Jonas Ries, EMBL, 2015'};
end