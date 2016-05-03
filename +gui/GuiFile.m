classdef GuiFile< interfaces.GuiModuleInterface & interfaces.LocDataInterface
    properties
        autosavetimer
        loaders
        savers
    end
    methods
        function obj=GuiFile(varargin)
            obj@interfaces.GuiModuleInterface(varargin{:})
            obj.outputParameters={'group_dx','group_dt'};
            obj.excludeFromSave={'filelist_long','loadmodule','savemodule'};
        end
        function initGui(obj)
            
            fw=obj.guiPar.FieldWidth;
            fh=obj.guiPar.FieldHeight;
            allloaders=pluginnames('File','Load');
            
            lp={};
            for k=1:length(allloaders)

                loaderhandle=uipanel(obj.handle,'Units','pixels','Position',[1,5.25*fh,2.9*fw,3.*fh],'Visible','off');
                obj.guihandles.(['loaderpanel_' num2str(k)])=loaderhandle;
                loaderp=plugin('File','Load',allloaders{k},loaderhandle,obj.P);
%                 infom.pluginpath={'File','Load',obj.loaders{k}};
                info=loaderp.info;
                lp{k}=info.name;
                loaderp.attachLocData(obj.locData);
                loaderp.makeGui;
                obj.children.(['loader_' num2str(k)])=loaderp;
                loaders{k}=loaderp;
            end
            obj.loaders=loaders;
            obj.loaders{1}.handle.Visible='on';
            obj.guihandles.loadmodule.String=lp;
            
            allsavers=pluginnames('File','Save');
            

%             obj.savers= plugintemp.plugins('File','Save',allsavers{1},[],obj.P); %to initialize
            for k=1:length(allsavers)
                saverhandle=uipanel(obj.handle,'Units','pixels','Position',[1,1*fh,2.9*fw,3.*fh],'Visible','off');
                obj.guihandles.(['saverpanel_' num2str(k)])=saverhandle;
                
                saver=plugin('File','Save',allsavers{k},saverhandle,obj.P);
%                 saver.pluginpath={'File','Save',allsavers{k}};
                saver.attachLocData(obj.locData);
                saver.makeGui;
%                 info=plugintemp.plugins('File','Save',obj.savers{k}).info;
                ls{k}=saver.info.name;
                obj.children.(['saver_' num2str(k)])=saver;
                savers{k}=saver;
            end 
            obj.savers=savers;
            obj.savers{1}.handle.Visible='on';
            obj.guihandles.savemodule.String=ls;
            obj.guihandles.savemodule.Value=1;
            obj.addSynchronization('filelist_long',obj.guihandles.filelist_long,'String');
            
            %make file menu
            makefilemenu(obj);
            
            f=getParentFigure(obj.handle);
            c=uicontextmenu(f);
            obj.guihandles.filelist_long.UIContextMenu=c;
            m1 = uimenu(c,'Label','info','Callback',{@menu_callback,obj});
            m1 = uimenu(c,'Label','remove','Callback',{@menu_callback,obj});
            
        end
     
        function loadbutton_callback(obj, handle,actiondata,isadd)
            %load data 
            p=obj.getAllParameters;
            fm=p.mainfile;   
            if isempty(fm)
                fm=p.filelist_long.selection;
            end
            path=fileparts(fm);          
            loader=obj.loaders{p.loadmodule.Value};
%             p=copyfields(p,loader.getGuiParameters);
%             loader=plugin('File','Load',obj.loaders{p.loadmodule.Value},[],obj.P);
%             loader.attachLocData(obj.locData);
%             loader.makeGui;
            try
                 ext=loader.info.extensions;
                 title=loader.info.dialogtitle;
            catch
                ext='*.*';
                title='format not specified'
            end
            [f,pfad]=uigetfile(ext,title,path,'MultiSelect','on');
            if pfad %file selected
                if ~iscell(f)
                    f={f};
                end
                
                loader.clear([pfad f{1}],isadd)
%                 [~,~,ext]=fileparts(f{1});

%                 [mode, emptylocs]=getfilemode([pfad f{1}]);
%                 if isadd || ~emptylocs
%                     obj.locData.empty('filter');
%                 else
%                     obj.locData.empty;
%                 endgl
                
                for k=1:length(f)
                    [~,~,ext]=fileparts(f{k});
                    if ~isempty(strfind(f{k},'_sml'))
                        obj.setPar('mainfile',[pfad filesep f{k}]);
                    end
                    obj.status(['load: ' f{k}])
                    drawnow
                    loadfiles(obj,loader,f{k},pfad)
                end

                obj.status('file loaded')

                initGuiAfterLoad(obj)
                autosavecheck_callback(0,0,obj)
            end 
        end
        
       
        
        function remove_callback(obj,callobj,handle,actiondata)
            removefile=get(obj.guihandles.filelist_long,'Value');
            floc=obj.locData.loc.filenumber;
            removeind=floc==removefile;
            moveup=floc>removefile;
            
            obj.locData.loc.filenumber(moveup)= obj.locData.loc.filenumber(moveup)-1;
            obj.locData.removelocs(removeind);
            obj.locData.files.filenumberEnd=obj.locData.files.filenumberEnd-1;
            obj.locData.files.file(removefile)=[];
            obj.locData.SE.removeFile(removefile);
            
            fl=obj.getPar('filelist_long').String;
            fs=obj.getPar('filelist_short').String;
            fl(removefile)=[];
            fs(removefile)=[];
            obj.locData.regroup;
            obj.locData.filter;
            obj.setPar('filelist_long',fl,'String');
            obj.setPar('filelist_short',fs,'String');
        end     
        
        function setGuiParameters(obj,p,varargin)
            setGuiParameters@interfaces.GuiModuleInterface(obj,p,varargin{:});
            if isfield(p,'filelist_long')
            obj.guihandles.filelist_long.String=p.filelist_long.String;   
            end
           
        end
        function group_callback(obj,a,b)
            obj.locData.regroup;
%             group_callback(0, 0,obj);
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function delete(obj)
            delete(obj.autosavetimer)
        end
        
    end
end

function savemenu_callback(object,event,obj,what)
switch what
    case 'sml'
        savemodule=1;
    case 'tif'
        savemodule=2;
end
obj.guihandles.savemodule.Value=savemodule;
save_callback(0,0,obj)
end

function save_callback(object,event,obj)
p=obj.getAllParameters;
saver=obj.savers{p.savemodule.Value};
psave=obj.getAllParameters(saver.inputParameters);
saver.save(psave);
end


function savemode_callback(data,b,obj)
for k=1:length(obj.savers)
    obj.savers{k}.handle.Visible='off';
end
obj.savers{data.Value}.handle.Visible='on';
end

function loadmode_callback(data,b,obj)
for k=1:length(obj.loaders)
    obj.loaders{k}.handle.Visible='off';
end
obj.loaders{data.Value}.handle.Visible='on';
end



function autosavecheck_callback(a,b,obj)
p=obj.getAllParameters;
t=obj.autosavetimer;
%creaste timer if empty or invalid
if isempty(t)||~isa(t,'timer')||~isvalid(t)
    t=timer;
    t.Period=p.autosavetime*60;
    t.StartDelay=t.Period;
    t.TimerFcn={@autosave_timer,obj};
    t.ExecutionMode='fixedRate';
    obj.autosavetimer=t;
end

if p.autosavecheck %deactivate
    if strcmpi(t.Running,'off')
    start(t);
    end
else
    if strcmpi(t.Running,'on')
    stop(t);
    end
end
end

function autosave_timer(a,b,obj)
p.mainGui=obj.getPar('mainGui');
p.saveroi=false;
if obj.guihandles.autosavecheck.Value
    savesml(obj.locData,'settings/temp/autosave_sml',p)
    time=datetime('now');
    disp(['autosave: ' num2str(time.Hour) ':' num2str(time.Minute)])
end
end

function autosavetime_callback(a,b,obj)
p=obj.getGuiParameters;
t=obj.autosavetimer;
if ~isempty(t)||isa(t,'timer')
    if strcmpi(t.Running,'on')
        stop(t);
    end
    obj.autosavetimer.Period=p.autosavetime*60;
    obj.autosavetimer.StartDelay=obj.autosavetimer.Period;
    if obj.guihandles.autosavecheck.Value
        start(t);
    end
end
end

function loadfiles(obj,loader,f,pfad)
mode=getfilemode([pfad f]);
if strcmp(mode,'tif')
    si=checkforsingleimages([pfad f]);             
    if si==1
        obj.setPar('filelist_localize',[pfad f]) %communication with localizer
        maing=obj.getPar('mainGui');
        maing.setmaintab(2);
        obj.locData.clear;
        return
    elseif si==0
%         return
    end
end


par=obj.getAllParameters(loader.inputParameters);
par=copyfields(par,loader.getGuiParameters);
loader.load(par,[pfad f]);
end

function menu_callback(menuobj,b,obj)
switch menuobj.Label
    case 'info'
        fin=obj.getPar('filelist_long').Value;
        info=obj.locData.files.file(fin).info;
        texta={};
        fn=fieldnames(info);
        for k=1:length(fn)
            ph=info.(fn{k});
            
            if iscell(ph)
                v='cell';
            elseif ischar(ph)
                v=ph;
            elseif isnumeric(ph)
                v=num2str(ph);
            end
                
            texta{end+1}=[fn{k} ': ' sprintf('\t') v];
        end
            
            listdlg('ListString',texta,'ListSize',[800,800]);
    case 'remove'
        obj.remove_callback;

end
end

function makefilemenu(obj)
hsmap=obj.getPar('mainGuihandle');
hfile=uimenu(hsmap,'Label','File');
hload=uimenu(hfile,'Label','load','Callback',{@obj.loadbutton_callback,0});
hsavesml=uimenu(hfile,'Label','save SML','Callback',{@savemenu_callback,obj,'sml'});
hsavetif=uimenu(hfile,'Label','save Tif','Callback',{@savemenu_callback,obj,'tif'});

tobj=findobj(hsmap.Children,'Label','SMAP');
itarget=find(hsmap.Children==tobj);
ithis=find(hsmap.Children==hfile);
indold=1:length(hsmap.Children);
indnew=[ setdiff(indold,[ithis itarget])  ithis itarget];

% hdummy=hsmap.Children(itarget);
% hsmap.Children([itarget(1) ithis(1)])=hsmap.Children([ithis(1) itarget(1)]);
hsmap.Children=hsmap.Children(indnew);
% hsmap.Children(ithis)=hdummy;
end

function pard=guidef(obj)
pard.load.object=struct('Style','pushbutton','String','Load','Callback',{{@obj.loadbutton_callback,0}});
pard.load.position=[4.75,2.5];
pard.load.Width=0.7;
pard.load.Height=1.5;

pard.add.object=struct('Style','pushbutton','String','Add','Callback',{{@obj.loadbutton_callback,1}});
pard.add.position=[4.75,3.2];
pard.add.Width=0.7;
pard.add.Height=1.5;

pard.loadmodule.object=struct('Style','popupmenu','String',{{'auto'}},'Callback',{{@loadmode_callback,obj}});
pard.loadmodule.position=[4.5,1];
 pard.loadmodule.Width=1.5;
 pard.loadmodule.object.TooltipString='select saver plugin';
 
% pard.updateGuiPar.object=struct('Style','checkbox','String','load Gui Parameters');
% pard.updateGuiPar.position=[5.5,1];
%  pard.updateGuiPar.Width=1.5;
% pard.updateGuiPar.TooltipString='Restore Gui parameters saved with localization data';

% 
% pard.remove.object=struct('Style','pushbutton','String','remove','Callback',{{@obj.remove_callback,1}});
% pard.remove.position=[4,4.5];
% pard.remove.Width=0.5;
% pard.add.Height=1.5;


pard.filelist_long.object=struct('Style','Listbox','String',{'x://'});
pard.filelist_long.position=[3,1];
pard.filelist_long.Width=4;
pard.filelist_long.Height=3;

pard.autosavecheck.object=struct('Style','checkbox','String','Auto save','Value',1);
pard.autosavecheck.position=[9,4];
pard.autosavecheck.Width=1;
pard.as_tdx.object=struct('Style','text','String','Interval (min)');
pard.as_tdx.position=[10,4];
pard.as_tdx.Width=0.65;
pard.autosavetime.object=struct('Style','edit','String','30','Callback',{{@autosavetime_callback,obj}});
pard.autosavetime.position=[10,4.65];
pard.autosavetime.Width=0.35;

pard.savemodule.object=struct('Style','popupmenu','String',{{'_sml','final image','raw images','_fitpos','settings'}},...
    'Callback',{{@savemode_callback,obj}});
pard.savemodule.position=[9,1.];
pard.savemodule.Width=1.5;

pard.save.object=struct('Style','pushbutton','String','Save','Callback',{{@save_callback,obj}});
pard.save.position=[9,2.5];
pard.save.Width=1;

pard.group_b.object=struct('Style','pushbutton','String','Group','Callback',{{@obj.group_callback}});
pard.group_b.position=[5,4];
pard.group_b.Width=1;

pard.group_tdx.object=struct('Style','text','String','dX (nm)');
pard.group_tdx.position=[6,4];
pard.group_tdx.Width=0.65;
pard.group_dx.object=struct('Style','edit','String','75');
pard.group_dx.position=[6,4.65];
pard.group_dx.Width=0.35;

pard.group_tdt.object=struct('Style','text','String','dT (frames)');
pard.group_tdt.position=[7,4];
pard.group_tdt.Width=0.65;
pard.group_dt.object=struct('Style','edit','String','1');
pard.group_dt.position=[7,4.65];
pard.group_dt.Width=0.35;

pard.syncParameters={{'filelist_long','filelist_long',{'String'}}};
pard.outputParameters= {'group_dx','group_dt'};
pard.inputParameters={'mainfile'};


pard.load.object.TooltipString='load localization data or image. Load at least one localization data before adding an image.';
pard.add.object.TooltipString='add localization data or image';
% pard.remove.object.TooltipString='remove file';
pard.savemodule.object.TooltipString='select what to save';
pard.group_b.object.TooltipString='group localizations in consecutive frames';
pard.group_dx.object.TooltipString=sprintf('distance in nm which two locs can be apart \n and still grouped together');
pard.group_dt.object.TooltipString=sprintf('number of frames locs can be missing \n and still grouped together');

end