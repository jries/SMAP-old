classdef Batchprocessor<interfaces.GuiModuleInterface&interfaces.LocDataInterface
    properties %(Access=private)
        mainbatchfile
    end
    properties

    end
    methods
        function obj=Batchprocessor(varargin)
            obj@interfaces.GuiModuleInterface(varargin{:})
            if isempty(obj.handle)||~isvalid(obj.handle)
                obj.handle=figure('Units','normalized','Units','pixels','Position',[150,200,900,300],'Name','Batch','NumberTitle','off');
                delete(obj.handle.Children);
            end
        end
        function initGui(obj)
            if ~isempty(obj.mainbatchfile)&&exist(obj.mainbatchfile,'file')
                obj.guihandles.mainbatchfile.String=obj.mainbatchfile;
                obj.guihandles.filelist.String={obj.mainbatchfile};
            end
        end
        function pard=pardef(obj)
            pard=pardef(obj);
        end
        function mainbatchfileb_callback(obj,a,b)
            p=obj.getGuiParameters;
            [f, path]=uigetfile(p.mainbatchfile);
            if f
                obj.guihandles.mainbatchfile.String=[path f];
                obj.mainbatchfile=[path f];
            end
        end
        function addb_callback(obj,a,b)
            p=obj.getGuiParameters;
            if ~isempty(p.filelist.selection)
                path=fileparts(p.filelist.selection);
            else
                path=fileparts(p.mainbatchfile);
            end
            [f,path]=uigetfile(fullfile(path,'*_batch.mat;*.tif'),'Select image or batch files');
            if ~f
                return
            end
            str=p.filelist.String;
            str{end+1}=[path f];
            obj.guihandles.filelist.String=str;
            
            if ~isempty(strfind(f,'_batch.mat'))&&~exist(p.mainbatchfile,'file') %put as main batch file
                obj.guihandles.mainbatchfile.String=[path f];
            end
            
        end
        function adddirb_callback(obj,a,b)
            p=obj.getGuiParameters;
            if ~isempty(p.filelist.selection)
                path=fileparts(p.filelist.selection);
            else
                path=fileparts(p.mainbatchfile);
            end
            [path]=uigetdirs(path,'Select image or batch directories');
            if isempty(path)
                return
            end
            
            str=p.filelist.String;
            for k=1:length(path)
                p.hstatus=obj.guihandles.status;
                imf=findimageindir(path{k},p);
                if ~isempty(imf)
                    for l=1:length(imf)
                        str{end+1}=imf{l};
                    end
                end
                    
            end
            obj.guihandles.filelist.String=str;          
        end
        function removeb_callback(obj,a,b)
            p=obj.getGuiParameters;
            vo=p.filelist.Value;
            str=p.filelist.String;
            str(vo)=[];
            obj.guihandles.filelist.String=str;
            obj.guihandles.filelist.Value=min(vo,length(str));
        end
        function processb_callback(obj,a,b)
            obj.setPar('loc_preview',false);
            p=obj.getGuiParameters;
            
            filelist=p.filelist.String;
            for k=1:length(filelist)
                filen=filelist{k};
                if ~isempty(strfind(filen,'.tif'))
                    obj.processtiff(filen);
                elseif ~isempty(strfind(filen,'.mat')) 
                    if p.useforall
                        obj.processtiffromWF(filen);
                    else
                        obj.processWF(filen);
                    end
                end
            end
        end
        
        function processWF(obj,wffile)
            wf=interfaces.Workflow([],obj.P);
            wf.attachLocData(obj.locData);
            wf.makeGui;
            wf.load(wffile);
            %workflow: recover tiff file and start
            tiffile=wf.module(wf.startmodule).getGuiParameters.tiffile;
            wf.module(wf.startmodule).addFile(tiffile);
            wf.run;   
        end
        function processtiff(obj,tiffile)
            p=obj.getGuiParameters;
            wf=interfaces.Workflow([],obj.P);
            wf.attachLocData(obj.locData);
            wf.makeGui;
            wf.load(p.mainbatchfile);
            wf.module(wf.startmodule).addFile(tiffile);
            wf.run;
        end
        
        function processtiffromWF(obj,wffile)
            wf=interfaces.Workflow([],obj.P);
            wf.attachLocData(obj.locData);
            wf.makeGui;
            wf.load(wffile);
            %workflow: recover tiff file and start
            tiffile=wf.module(wf.startmodule).getGuiParameters.tiffile;
            obj.processtiff(tiffile);
        end
        
    end

end

function img=findimageindir(path,p)
img={};
files=dir([path filesep '*.tif']);
if ~isempty(files)
    img={[path filesep files(1).name]};
    return
end
files=dir([path filesep 'Pos*']);
for k=1:length(files)
    if files(k).isdir
        path2=[path filesep files(k).name filesep];
        files2=dir([path2 '*.tif']);
        if ~isempty(files2)
            img={[path2 files2(1).name]};
            return
        end
    end
end
searchstr=p.adddir_mask;
mintiffs=p.adddir_minimages;
files=dir([path filesep searchstr]);
for k=1:length(files)
    if files(k).isdir 
        tiffiles=(dir([path filesep files(k).name filesep '*.tif']));
        addpathh='';
        if isempty(tiffiles)
            tiffiles=(dir([path filesep files(k).name filesep 'Pos0' filesep '*.tif']));
            addpathh=['Pos0' filesep];
        end
        if length(tiffiles)>mintiffs
            img{end+1}=[path filesep files(k).name filesep addpathh tiffiles(1).name];
        end
    end 
    p.hstatus.String=['directory ' num2str(k) ' of ' num2str(length(files))];
    drawnow;
end
end

function pard=pardef(obj)
pard.filelist.object=struct('Style','listbox');%,'Callback',@obj.moduleselect_callback);
pard.filelist.position=[10,1];
pard.filelist.Height=8;
pard.filelist.Width=3;

pard.mainbatchfile.object=struct('Style','edit','String','x://*_batch.mat');
pard.mainbatchfile.position=[1,1];
pard.mainbatchfile.Width=3;

pard.mainbatchfile_button.object=struct('Style','pushbutton','String','load main batchfile','Callback',@obj.mainbatchfileb_callback);
pard.mainbatchfile_button.position=[1,4];
pard.mainbatchfile_button.Width=1;

pard.add_button.object=struct('Style','pushbutton','String','add','Callback',@obj.addb_callback);
pard.add_button.position=[3,4];
pard.add_button.Width=1;

pard.adddir_button.object=struct('Style','pushbutton','String','add directories','Callback',@obj.adddirb_callback);
pard.adddir_button.position=[4,4];
pard.adddir_button.Width=1;

pard.adddir_mask.object=struct('Style','edit','String','*_time_*');
pard.adddir_mask.position=[5,4];
pard.adddir_mask.Width=1;

pard.adddir_t.object=struct('Style','text','String','> #images');
pard.adddir_t.position=[6,4];
pard.adddir_t.Width=.5;
pard.adddir_minimages.object=struct('Style','edit','String','10');
pard.adddir_minimages.position=[6,4.5];
pard.adddir_minimages.Width=.5;

pard.remove_button.object=struct('Style','pushbutton','String','remove','Callback',@obj.removeb_callback);
pard.remove_button.position=[8,4];
pard.remove_button.Width=1;

pard.process_button.object=struct('Style','pushbutton','String','Batch process','Callback',@obj.processb_callback);
pard.process_button.position=[10,4];
pard.process_button.Height=2;

pard.useforall.object=struct('Style','checkbox','String','use for all','Value',0);
pard.useforall.position=[2,4];
pard.useforall.Height=1;

pard.status.object=struct('Style','text','String','status','Value',0);
pard.status.position=[12,1];
pard.status.Width=4;
end