classdef driftcorrectionXYZ<interfaces.DialogProcessor
    methods
        function obj=driftcorrectionXYZ(varargin)        
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.inputParameters={'layer1_'};
                obj.history=true;
                obj.guiselector.show=true;
        end
        
        function out=run(obj,p)
            
            %batch: many sml files loaded:
                %do it per file, save only file.
                %locData.copy, remove other filenames from .loc, .grouploc,
                %. files.file(k)
                %save copy
                obj.setPar('undoModule','driftfeature');
            notify(obj.P,'backup4undo');
            groupcheck=obj.locData.isgrouped(1);
            numberOfFiles=obj.locData.files.filenumberEnd;
            if p.drift_individual&&numberOfFiles>1
                for k=1:numberOfFiles
                    lochere=obj.locData.copy;
                    lochere.files.fileNumberEnd=1;
                    lochere.files.file=lochere.files.file(k);
                    badind=lochere.loc.filenumber~=k;
                    lochere.removelocs(badind);
                    lochere.regroup;
                    lochere.loc.filenumber=lochere.loc.filenumber*0+1;
                    
                    locs=lochere.getloc({'frame','xnm','ynm'},'position','all','grouping',groupcheck);
                    p.maxframeall=max(lochere.loc.frame);
                    p.framestart=1;
                    p.framestop=p.maxframeall;
                    [drift,driftinfo]=finddriftfeature(locs,p);
                    locsall=copyfields([],lochere.loc,{'xnm','ynm','frame','filenumber'});
                    locsnew=applydriftcorrection(drift,locsall);
                    lochere.loc=copyfields(lochere.loc,locsnew,{'xnm','ynm'});
                    lochere.files.file(1).driftinfo=driftinfo;
                    obj.locData.files.file(k).driftinfo=driftinfo;
                    fn=lochere.files.file(1).name;
                    fnn=strrep(fn,'_sml','_driftc_sml');
                    lochere.savelocs(fnn); 
                    obj.locData.loc.xnm(~badind)=lochere.loc.xnm;
                    obj.locData.loc.ynm(~badind)=lochere.loc.ynm;
                end
                obj.locData.regroup;
            else
                locs=obj.locData.getloc({'frame','xnm','ynm','znm'},'position','fov','grouping',groupcheck);
                if length(locs.xnm)<100
                    locs=obj.locData.getloc({'frame','xnm','ynm','znm'},'position','all');
                end
                p.maxframeall=max(obj.locData.loc.frame);
                p.framestart=p.layer1_.frame_min;
                p.framestop=min(p.layer1_.frame_max,p.maxframeall);
                [drift,driftinfo]=finddriftfeature(locs,p);
                if ~isempty(locs.znm)&&p.correctz
                    [driftz,driftinfoz]=finddriftfeatureZ(locs,p);
                    drift=copyfields(drift,driftz,'z');
                     driftinfo=copyfields(driftinfo,driftinfoz);
                     fieldc={'xnm','ynm','znm'};
                else
                    fieldc={'xnm','ynm'};
                end
                locsall=copyfields([],obj.locData.loc,{fieldc{:},'frame','filenumber'});
%                 grouplocsall=copyfields([],obj.locData.grouploc,{'xnm','ynm','frame','filenumber'});
                locsnew=applydriftcorrection(drift,locsall);
%                 grouplocsnew=applydriftcorrection(drift,grouplocsall);
                obj.locData.loc=copyfields(obj.locData.loc,locsnew,fieldc);
%                 obj.locData.grouploc=copyfields(obj.locData.grouploc,grouplocsnew,{'xnm','ynm'});

                if length(unique(locsall.filenumber))>1
                    locsnew.channel=locsall.filenumber;
%                     grouplocsnew.channel=grouplocsall.filenumber;
                    obj.locData.loc=copyfields(obj.locData.loc,locsnew,{'channel'});
%                     obj.locData.grouploc=copyfields(obj.locData.grouploc,grouplocsnew,{'channel'});
                end

               obj.locData.files(obj.locData.loc.filenumber(1)).file.driftinfo=driftinfo;
                fn=obj.locData.files(obj.locData.loc.filenumber(1)).file.name;
                fnn=strrep(fn,'_sml','_driftc_sml');
    %             fnn=[fn(1:end-7) '_driftc_sml.mat'];
                obj.locData.savelocs(fnn);   
                obj.locData.regroup;
                out=1;
            end
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef
pard.texta.object=struct('String','timepoints','Style','text');
pard.texta.position=[1,1];


pard.drift_timepoints.object=struct('String','10','Style','edit');
pard.drift_timepoints.object.TooltipString=sprintf('whole data is divided into timepoints individual \n blocks. Range: 7-20');
pard.drift_timepoints.position=[1,2];
pard.drift_timepoints.isnumeric=1;

pard.text1.object=struct('String','pixrec nm','Style','text');
pard.text1.position=[2,1];
pard.text1.Optional=true;

pard.drift_pixrec.object=struct('String','15','Style','edit');
pard.drift_pixrec.position=[2,2];
% pard.drift_pixrec.isnumeric=1;
pard.drift_pixrec.object.TooltipString=sprintf('pixel size (nm) for reconstruction. \n Smaller for well defined peak. But slower \n Range: 10-25');
pard.drift_pixrec.Optional=true;

pard.text2.object=struct('String','window pix','Style','text');
pard.text2.position=[3,1];
pard.text2.Optional=true;

pard.drift_window.object=struct('String','11','Style','edit');
pard.drift_window.position=[3,2];
% pard.drift_window.isnumeric=1;
pard.drift_window.object.TooltipString=sprintf('size of region for peakfinding (ellipt. Gaussian). \n should be small, but cover clear maximum. \n Range: 7-15');
pard.drift_window.Optional=true;

pard.text3.object=struct('String','maxdrift nm','Style','text');
pard.text3.position=[4,1];
pard.text3.Optional=true;

pard.drift_maxdrift.object=struct('String','1000','Style','edit');
pard.drift_maxdrift.position=[4,2];

pard.drift_maxdrift.object.TooltipString=sprintf('Maximum drift expected. \n Smaller if data is sparse and wrong peak found. \n larger if no clear peak found. \n Range 250-2000');
pard.drift_maxdrift.Optional=true;

pard.correctz.object=struct('String','Correct z-drift','Style','checkbox','Value',0);
pard.correctz.position=[1,3];


pard.drift_pixreczt.object=struct('String','z binwidth nm','Style','text');
pard.drift_pixreczt.position=[2,3];
pard.drift_pixreczt.Optional=true;

pard.drift_pixrecz.object=struct('String','5','Style','edit');
pard.drift_pixrecz.position=[2,4];
pard.drift_pixrecz.isnumeric=1;
pard.drift_pixrecz.object.TooltipString=sprintf('pixel size (nm) for reconstruction. \n Smaller for well defined peak. But slower \n Range: 10-25');
pard.drift_pixrecz.Optional=true;

pard.drift_windowzt.object=struct('String','z fit window pix','Style','text');
pard.drift_windowzt.position=[3,3];
pard.drift_windowzt.Optional=true;

pard.drift_windowz.object=struct('String','15','Style','edit');
pard.drift_windowz.position=[3,4];
pard.drift_windowz.isnumeric=1;
pard.drift_windowz.object.TooltipString=sprintf('size of region for peakfinding (ellipt. Gaussian). \n should be small, but cover clear maximum. \n Range: 7-15');
pard.drift_windowz.Optional=true;

pard.zranget.object=struct('String','zrange nm','Style','text');
pard.zranget.position=[4,3];
pard.zranget.Optional=true;

pard.zrange.object=struct('String','-400 400','Style','edit');
pard.zrange.position=[4,4];
pard.zrange.Optional=true;

pard.slicewidtht.object=struct('String','slice width nm','Style','text');
pard.slicewidtht.position=[5,3];
pard.slicewidtht.Optional=true;

pard.slicewidth.object=struct('String','200','Style','edit');
pard.slicewidth.position=[5,4];
pard.slicewidth.Optional=true;



pard.drift_reference.object=struct('String','reference is last frame','Style','checkbox');
pard.drift_reference.position=[6,1];
pard.drift_reference.Optional=true;
pard.drift_reference.object.TooltipString=sprintf('If checked, drift at end of data set is set to zero. \n Useful for sequential acquisition, use this for first data set.');
pard.drift_reference.Width=2;
pard.drift_reference.Optional=true;

pard.drift_individual.object=struct('String','correct every file individually','Style','checkbox','Value',1);
pard.drift_individual.position=[7,1];
pard.drift_individual.Width=2;
pard.drift_individual.Optional=true;

pard.plugininfo.name='drift correctiom X,Y,Z';
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description={'Drift correction based on cross-correlation.','Algorithm: the data set is divided into [timepoints] blocks, for which superresolution images are calculated. The displacement between all images is calcualted with a FFT-based cross-correlation algorithm. The position of the maxima of the cross-correlation curve are fitted with sub-pixel accuracy with a free elliptical Gaussian.',...
    'A robust estimator is used to calculate the drift vs frame from all pairwise displacements.','All localiaztions visible in the superresolution image are used to infer the drift. Use [Render]...[Layer] to control this.',...
    'If two files are loaded, their drift is calculated together and they are saved as one file with their filenumbers copied to the channel field.',' ','(c) Jonas Ries, EMBL, 2015'};
end