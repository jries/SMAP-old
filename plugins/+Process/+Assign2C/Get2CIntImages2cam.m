classdef Get2CIntImages2cam<interfaces.DialogProcessor
    % gets intensities from camera images at positions of localizations and
    % at transformed positions
    properties (Access=private)
        figure
    end
    methods
        function obj=Get2CIntImages2cam(varargin)   
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.history=true;
            obj.showresults=false;
        end
        
        function out=run(obj,p)
            out=[];
            if isempty(obj.figure)||~isvalid(obj.figure)
                obj.figure=figure;
            end
            f=obj.figure;
            f.Visible='on';
            wffile='settings/workflows/get2CIntensityImagesWF_group.mat';
            wf=interfaces.Workflow(f,obj.P);
            wf.attachLocData(obj.locData);
            wf.makeGui;
            wf.load(wffile);
            
      
            p.loaderclass=wf.module('TifLoader').getSingleGuiParameter('loaderclass'); 
            p.loaderclass.Value=1;%single MM
           
             wf.module('TifLoader').setGuiParameters(p);
                  
               p.framestop=max(obj.locData.loc.frame);

            file=obj.locData.files.file;
            wf.module('IntLoc2posN').filestruc=file;
            wf.module('IntLoc2posN').setGuiParameters(p);
            
            p.loc_ROIsize=11;
            p.loc.fitgrouped=true;
            wf.module('RoiCutterWF_groupExt').setGuiParameters(p);

            pe=obj.children.evaluate.getGuiParameters(true);
            wf.module('EvaluateIntensity_s').setGuiParameters(pe,true);
            obj.setPar('loc_preview',false);


            
            % now first to ref, then do target. Later: if files are same:
            % do at the same time to save time...
            if p.evaltarget
                wf.module('TifLoader').addFile(p.tiffiletarget,true);   
                wf.module('EvaluateIntensity_s').extension='t';
                wf.module('IntLoc2posN').setGuiParameters(struct('transformtotarget',true));
                wf.run;
            end
            if p.evalref
                wf.module('TifLoader').addFile(p.tiffileref,true);   
                wf.module('EvaluateIntensity_s').extension='r';
                
                wf.module('IntLoc2posN').setGuiParameters(struct('transformtotarget',false));
                wf.run;
            end
            delete(f);
                        fo=strrep(obj.locData.files.file(1).name,'_sml.mat','_dc_sml.mat');
            obj.locData.savelocs(fo);
            obj.locData.regroup;
            obj.setPar('locFields',fieldnames(obj.locData.loc))

        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function initGui(obj)
            par.Vpos=3;
            par.Xpos=3;
            obj.children.evaluate=makeplugin(obj,{'WorkflowModules','IntensityCalculator','EvaluateIntensity'},obj.handle,par);
%             obj.guihandles.loadbutton.Callback=@obj.loadbutton;

        end
        function loadbutton_T(obj,a,b)
            fn=obj.guihandles.Tfile.String;
            [f,path]=uigetfile(fn,'Select transformation file _T.mat');
            if f
                obj.guihandles.Tfile.String=[path f];
            end      
        end
        function loadbutton_tif(obj,a,b,field)
            fn= obj.guihandles.(field).String;
            if isempty(fn)
                fn=getrawtifpath(obj.locData);
            end
            
            [f,path]=uigetfile(fn,'Select raw tiff file');
            if f
                obj.guihandles.(field).String=[path f];
            end      
        end
    end
end




function pard=guidef(obj)
% 
pard.evaltarget.object=struct('Style','checkbox','String','target','Value',1);
pard.evaltarget.position=[1,1];
pard.evaltarget.Width=0.7;

pard.tiffiletarget.object=struct('Style','edit','String','');
pard.tiffiletarget.position=[1,1.7];
pard.tiffiletarget.Width=2.6;

pard.loadbuttontiftarget.object=struct('Style','pushbutton','String','load tif','Callback',{{@obj.loadbutton_tif,'tiffiletarget'}});
pard.loadbuttontiftarget.position=[1,4.3];
pard.loadbuttontiftarget.Width=0.7;

pard.evalref.object=struct('Style','checkbox','String','reference');
pard.evalref.position=[2,1];
pard.evalref.Width=0.7;

pard.tiffileref.object=struct('Style','edit','String','');
pard.tiffileref.position=[2,1.7];
pard.tiffileref.Width=2.6;

pard.loadbuttontifref.object=struct('Style','pushbutton','String','load tif','Callback',{{@obj.loadbutton_tif,'tiffileref'}});
pard.loadbuttontifref.position=[2,4.3];
pard.loadbuttontifref.Width=0.7;
% pard.Tmode.object=struct('Style','popupmenu','String',{{'target','reference'}});
% pard.Tmode.position=[2,1];
% pard.Tmode.Width=.7;
% pard.Tmode.TooltipString=sprintf('transformation file. Created e.g. with RegisterLocs');

pard.Tfile.object=struct('Style','edit','String','');
pard.Tfile.position=[3,1.7];
pard.Tfile.Width=2.6;
pard.Tfile.TooltipString=sprintf('transformation file. Created e.g. with RegisterLocs');

pard.loadbuttonT.object=struct('Style','pushbutton','String','load T','Callback',@obj.loadbutton_T);
pard.loadbuttonT.position=[3,4.3];
pard.loadbuttonT.Width=0.7;
pard.loadbuttonT.TooltipString=pard.Tfile.TooltipString;

pard.mirrorem.object=struct('Style','checkbox','String','EM mirror','Value',1);
pard.mirrorem.position=[4,1];
pard.mirrorem.Optional=true;
pard.mirrorem.Width=0.8;

pard.syncParameters={{'transformationfile','Tfile',{'String'}}};
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description=sprintf(['This plugin gets intensities from camera images at positions of localizations and at transformed positions \n',...
    'This plugin uses a transformation to find for every localization the position in the other channel and then determines the intensity in both channels.\n',...
    '1.	Load a transformation\n',...
    '2.	Per default, this plugin does median filtering. Select the spatial and temporal spacing for this (dx, dt).\n',...
    '3.	Select one or several plugins which determine the intensity:\n',...
    '\t a.	Roi2int_sum: uses a ROI (set size) to determine intensity, and a larger ROI for the background.\n',...
    '\t b.	Roi2int_fit: Uses a Gaussian fit to determine intensity and background. The position is fixed to the fitted position. You can use the fitted PSF size or fix it. If fit on BG is checked, the background is subtracted prior to fitting and the fit is performed with background set to zero. Otherwise the background is a fitting parameter.\n',...
    '4.	Press Run and when asked select the original raw camera images. The results are automatically saved with the _dc in the file name.\n']);
pard.plugininfo.name='2C intensities from images 2 cam';
% pard.plugininfo.description
end



function tiffile=getrawtifpath(locData)
    if isfield(locData.files.file(1).info,'imagefile')
        tiffile=locData.files.file(1).info.imagefile;
    else
        tiffile=locData.files.file(1).info.basefile;
    end
    if ~exist(tiffile,'file')
        disp('Get2CIntImagesWF ine 40: check if it works')
        tiffile=strrep(tiffile,'\','/');
        ind=strfind(tiffile,'/');
        for k=1:length(ind)
            tiffileh=[tiffile(1:ind(k)) '_b_' tiffile(ind(k)+1:end)];
            if exist(tiffileh,'file')
                tiffile=tiffileh;
            end
        end
    end
    if ~exist(tiffile,'file')
        tiffile=strrep(locData.files.file(1).name,'_sml.mat','.tif');
    end
end

function plugino=makeplugin(obj,pluginpath,handle,pgui)
plugino=plugin(pluginpath{:});
plugino.attachPar(obj.P);
plugino.setGuiAppearence(pgui);
plugino.attachLocData(obj.locData);
plugino.handle=handle;
plugino.makeGui;
pluginname=pluginpath{end};
obj.children.(pluginname)=plugino;
end