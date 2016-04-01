classdef GuiLocalize<interfaces.GuiModuleInterface&interfaces.LocDataInterface
    properties
        mainworkflow
        customworkflows
        batchfile
    end
    methods
        function obj=GuiLocalize(varargin)
            obj@interfaces.GuiModuleInterface(varargin{:})
        end
        function makeGui(obj)
            h.loctab=uitabgroup(obj.handle,'Position',[0 .08 1 .92]);

            h.frame=uitab(h.loctab,'Title','Input Image');
            h.filter=uitab(h.loctab,'Title','Peak Finder');
            h.fit=uitab(h.loctab,'Title','Fitter');
            h.locprocess=uitab(h.loctab,'Title','Localizations');
%             h.workflow=uitab(h.loctab,'Title','Workflow');
            
            h.previewbutton=uicontrol(obj.handle,'Style','pushbutton','String','Preview','Position',[10 2, 70 obj.guiPar.FieldHeight*1.3],...
                'FontSize',obj.guiPar.fontsize,'Callback',{@preview_callback,obj});
             h.previewframe=uicontrol(obj.handle,'Style','edit','String','1','Position',[180 2, 60 obj.guiPar.FieldHeight*1],...
                'FontSize',obj.guiPar.fontsize,'Callback',{@previewframe_callback,obj});
           h.previewframeslider=uicontrol(obj.handle,'Style','slider','Position',[80 2, 100 20],...
                'FontSize',obj.guiPar.fontsize,'Callback',{@previewframeslider_callback,obj});
            
            h.localizebutton=uicontrol(obj.handle,'Style','pushbutton','String','Localize','Position',[380 2, 100 obj.guiPar.FieldHeight*1.3],...
                'FontSize',obj.guiPar.fontsize,'Callback',{@localize_callback,obj});
            h.batchbutton=uicontrol(obj.handle,'Style','pushbutton','String','Batch','Position',[260 2, 70 obj.guiPar.FieldHeight*1.3],...
                'FontSize',obj.guiPar.fontsize,'Callback',{@batch_callback,obj});


            outputfig=figure(207);
            outputfig.Visible='off';
            obj.setPar('loc_outputfig',outputfig)

            tabsizeh=obj.guiPar.tabsize2;
            tabsizeh(4)=250;
            h.framepanel=uipanel(h.frame,'Unit','pixels','Position',tabsizeh); 
            h.filterpanel=uipanel(h.filter,'Unit','pixels','Position',tabsizeh);
            h.fitpanel=uipanel(h.fit,'Unit','pixels','Position',tabsizeh); 
            h.locpanel=uipanel(h.locprocess,'Unit','pixels','Position',tabsizeh); 
            obj.guihandles=h;
            
            replacestruct={{'hframe',h.framepanel},{'hfilter',h.filterpanel},{'hfit',h.fitpanel},{'hloc',h.locpanel}};
            obj.createGlobalSetting('mainLocalizeWFFile','Directories','Description file for fitting workflow, e.g. settings/fit_tif_wavelet.txt',struct('Style','file','String','settings/fit_tif_wavelet.txt'))
%             settingsfile='settings/mainSMLMLocalizeWF.txt';
            settingsfile=obj.getGlobalSetting('mainLocalizeWFFile');
            par=readstruct(settingsfile,replacestruct);
            if isempty(par)
                warndlg('cannot find settings file for fit workflow. Please set in menu SMAP/Preferences')
            end
            wffile=par.all.file;
            
%             [wffile,par]=mainSMLMLocalizeWF(h.framepanel,h.filterpanel,h.fitpanel);

            mainworkflow=interfaces.Workflow([],obj.P);
            mainworkflow.attachLocData(obj.locData);
            mainworkflow.setGuiAppearence(par.all)
            mainworkflow.processorgui=false;
            mainworkflow.makeGui;
            mainworkflow.load(wffile,par);
            obj.mainworkflow=mainworkflow; 
            obj.children.mainworkflow=mainworkflow;
            
            f=getParentFigure(obj.handle);
            c=uicontextmenu(f);
            h.loctab.UIContextMenu=c;
            m1 = uimenu(c,'Label','remove','Callback',{@menu_callback,obj});
            m3 = uimenu(c,'Label','add workflow','Callback',{@menu_callback,obj});
            
%             parwf.Vpos=3;
%             parwf.Xpos=1;
%             parwf.FieldHeight=obj.guiPar.FieldHeight-4;
%             workflow=interfaces.Workflow(h.workflow,obj.P);
%             workflow.attachLocData(obj.locData);
%             workflow.processorgui=false;
%             workflow.setGuiAppearence(parwf)
%             workflow.makeGui;
%             obj.customworkflow=workflow;
%             obj.children.workflow=workflow;
            
            
            previewframe_callback(0,0,obj)
            obj.addSynchronization('loc_fileinfo',[],[],@obj.update_slider)
            
            h.batchprocessor=uicontrol(h.framepanel,'Style','pushbutton','String','Batch Processor','Position',[0, 0, 150 obj.guiPar.FieldHeight*1],...
                'FontSize',obj.guiPar.fontsize,'Callback',{@batchprocessor_callback,obj});
            
        end
        
        function update_slider(obj,a,b)        
            fi=obj.getPar('loc_fileinfo');
            nf=fi.numberOfFrames;
            if isempty(nf)||isinf(nf)
                nf=obj.getPar('loc_previewframe')+1;
            end
            obj.guihandles.previewframeslider.Min=1;
            obj.guihandles.previewframeslider.Max=nf;
            pvf=min(nf,obj.getPar('loc_previewframe'));
            obj.guihandles.previewframeslider.Value=pvf;
            obj.setPar('loc_previewframe',pvf)
            obj.guihandles.previewframe.String=num2str(pvf);
%             obj.setPar('loc_previewframe',round(pf));
            obj.guihandles.previewframeslider.SliderStep=[ceil(nf/50) ceil(nf/50)*10]/nf;
        end
    end
end

function menu_callback(callobj,b,obj)
switch callobj.Label
    case 'remove'
        selected=obj.guihandles.loctab.SelectedTab;
        title=selected.Title;
        if strcmp(title(1:2),'WF')
            number=str2double(title(3:end));
            delete(obj.customworkflows{number});
            obj.customworkflows(number)=[];
            delete(selected);
        else
            disp('only custom workflows are deletable')
        end

    case 'add workflow'
        nwf=length(obj.customworkflows);
        name=['WF' num2str(nwf+1)];
        obj.guihandles.(['tab_' name])=uitab(obj.guihandles.loctab,'Title',name);
        tabsizeh=obj.guiPar.tabsize2;
        tabsizeh(4)=250;
        obj.guihandles.(['panel_' name])=uipanel(obj.guihandles.(['tab_' name]),'Unit','pixels','Position',tabsizeh); 
        module=interfaces.Workflow;
        module.processorgui=false;
        module.handle=obj.guihandles.(['panel_' name]);
        module.attachPar(obj.P);
        module.attachLocData(obj.locData);
        p.Vrim=3;
        p.Xrim=4;
        p.FieldHeight=obj.guiPar.FieldHeight-4;
        module.setGuiAppearence(p)
        module.makeGui;
        obj.children.(name)=module;
        if nwf==0
            obj.customworkflows={module};
        else
            obj.customworkflows{end+1}=module;
        end
        
        obj.guihandles.loctab.SelectedTab= obj.guihandles.(['tab_' name]);
end
end

function preview_callback(a,b,obj)
selected=(obj.guihandles.loctab.SelectedTab.Title);
if strcmp(selected(1:2),'WF')
    number=str2double(selected(3:end));
    startmodule=obj.customworkflows{number};
else
    startmodule=obj.mainworkflow;
end
obj.setPar('loc_preview',true)
% obj.globpar.parameters.preview=true;
% notify(obj.P,'loc_initialize')
% startmodule.initialize;
startmodule.run;
obj.status('preview done');
end

function localize_callback(a,b,obj)
selected=(obj.guihandles.loctab.SelectedTab.Title);
if strcmp(selected(1:2),'WF')
    number=str2double(selected(3:end));
    startmodule=obj.customworkflows{number};
else
    startmodule=obj.mainworkflow;
end
obj.setPar('loc_preview',false)
% obj.globpar.parameters.preview=false;
% notify(obj.P,'loc_initialize')
% startmodule.initialize;
startmodule.run;
obj.locData.regroup;
obj.setPar('locFields',fieldnames(obj.locData.loc));
obj.status('fitting done');
maingui=obj.getPar('mainGui');
maingui.setmaintab(3);
end

function batch_callback(a,b,obj)

[f,p]=uiputfile([fileparts(obj.getPar('filelist_localize')) '_batch.mat']);
if ~f
    return;
end

if strcmp(obj.guihandles.loctab.SelectedTab.Title,'Workflow')
    wf=obj.customworkflow;
else
    wf=obj.mainworkflow;
end

wf.save([p f])
obj.batchfile=[p f];
end

function previewframe_callback(a,b,obj)
pf=str2double(obj.guihandles.previewframe.String);
obj.setPar('loc_previewframe',round(pf));
if pf>1&&pf<obj.guihandles.previewframeslider.Max
obj.guihandles.previewframeslider.Value=round(pf);
end
end

function previewframeslider_callback(a,b,obj)
obj.guihandles.previewframe.String=num2str(round(obj.guihandles.previewframeslider.Value));
obj.setPar('loc_previewframe',round(obj.guihandles.previewframeslider.Value));
end

function plugino=makeplugin(obj,pluginpath,handle,pgui)
% p=obj.guiPar;
% p=copyfields(p,pgui);
plugino=plugin(pluginpath{:});
plugino.attachPar(obj.P);

% plugin.addGuiParameters(p);
plugino.attachLocData(obj.locData);

plugino.setGuiAppearence(pgui);
plugino.handle=handle;
plugino.makeGui;

pluginname=pluginpath{end};
obj.children.(pluginname)=plugino;
end

function batchprocessor_callback(a,b,obj)
batchprocessor=WorkflowModules.Batchprocessor([],obj.P);
if ~isempty(obj.batchfile)
    batchprocessor.mainbatchfile=obj.batchfile;
end
batchprocessor.attachLocData(obj.locData);
batchprocessor.makeGui;
end