classdef Get2CIntImagesWF<interfaces.DialogProcessor
    methods
        function obj=Get2CIntImagesWF(varargin)   
            obj@interfaces.DialogProcessor(varargin{:}) ;
        end
        
        function out=run(obj,p)
            
            f=figure(389);
            f.Visible='off';
            wffile='plugins/workflows/get2CIntensityImagesWF.mat';
            wf=interfaces.Workflow(f,obj.P);
            wf.attachLocData(obj.locData);
            wf.makeGui;
            wf.load(wffile);
            
            
            load(p.Tfile)
            file=obj.locData.files.file(p.dataselect.Value);
            fo=strrep(file.name,'_sml.mat','_dc_sml.mat');

            path=fileparts(file.name); %rather in top class, pass on
            [filen, path]=uigetfile([path filesep '*.tif'],file.name);
            tiffile=[path filen];
%             tiffile='/Users/ries/Documents/Data/3Ddc/MTActin/02_MT680_phalloidin647_1/img_000039971_Default_000.tif';
            wf.module('TifLoader').addFile(tiffile);

             
            p.loc_bg_dt=p.filtert;
            p.loc_bg_dx=p.filterx;
            p.loc_subtractbg=true;
            wf.module('MedianBGcalculator').setGuiParameters(p);

            wf.module('IntLoc2pos').filestruc=file;
            wf.module('IntLoc2pos').setGuiParameters(p);
%             obj.children.loc2pos.prerun
            pe=obj.children.evaluate.getGuiParameters(true);
            wf.module('EvaluateIntensity').setGuiParameters(pe,true);
            obj.setPar('loc_preview',false);
            wf.run;

            obj.locData.savelocs(fo);
            obj.locData.regroup;
            obj.setPar('locFields',fieldnames(obj.locData.loc))

        end
        function pard=pardef(obj)
            pard=pardef;
        end
        function initGui(obj)
            par.Vpos=4;
            par.Xpos=3;
            obj.children.evaluate=makeplugin(obj,{'WorkflowModules','IntensityCalculator','EvaluateIntensity'},obj.handle,par);
            obj.guihandles.loadbutton.Callback=@obj.loadbutton;

        end
        function loadbutton(obj,a,b)
            fn=obj.guihandles.Tfile.String;
            [f,path]=uigetfile(fn,'Select transformation file _T.mat');
            if f
                obj.guihandles.Tfile.String=[path f];
            end      
        end
    end
end




function pard=pardef

pard.dataselect.object=struct('Style','popupmenu','String','File');
pard.dataselect.position=[1,1];

%sum
% pard.checksum.object=struct('Style','checkbox','String','sum','Value',1);
% pard.checksum.position=[1,2];
pard.t1.object=struct('Style','text','String','Background filter');
pard.t1.position=[3,1];
pard.t1.Width=1;
pard.t2.object=struct('Style','text','String',' \xi');
pard.t2.position=[4,1];
pard.t2.Width=0.5;
pard.filterx.object=struct('Style','edit','String','3');
pard.filterx.position=[4,1.5];
pard.filterx.Width=0.5;

pard.t3.object=struct('Style','text','String','BG \tau');
pard.t3.position=[5,1];
pard.t3.Width=0.5;
pard.filtert.object=struct('Style','edit','String','100');
pard.filtert.position=[5,1.5];
pard.filtert.Width=0.5;

pard.Tfile.object=struct('Style','edit','String','settings/temp/temp_T.mat');
pard.Tfile.position=[8,1];
pard.Tfile.Width=3;

pard.loadbutton.object=struct('Style','pushbutton','String','load');
pard.loadbutton.position=[8,4];

pard.syncParameters={{'filelist_short','dataselect',{'String'}},{'transformationfile','Tfile',{'String'}}};
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