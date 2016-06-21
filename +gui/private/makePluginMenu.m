function [pout,hsmap]=makePluginMenu(obj,handle)
hsmap=uimenu(handle,'Label','SMAP');
hinfo=uimenu(hsmap,'Label','About SMAP...','Callback',@info_callback);
hglobalSettings=uimenu(hsmap,'Label','Preferences...','Callback',{@globalsettings_callback,obj});
hexit=uimenu(hsmap,'Label','Quit SMAP','Callback',{@exit_callback,obj});

hrename=uimenu(hsmap,'Label','Rename window','Callback',{@renamewindow_callback,obj});

hsimplegui=uimenu(hsmap,'Label','Hide advanced controls','Callback',{@simplegui_callback,obj});

hmainplugin=uimenu(handle,'Label','Plugins');
% hwf=uimenu(handle,'Label','Workflows');
mwf1=uimenu(hmainplugin,'Label','New workflow','Callback',{@makeplugin,obj,'Workflow'});
obj.loadGlobalSettings;
% names1={'File','Analyze','Process','Siteexplorer'};
names1=pluginnames;
nomenutypes={'WorkflowModule', 'WorkflowFitter','ROI_Evaluate','WorkflowIntensity'}; %dont put those processors into menu
for k=1:length(names1)
    names2=pluginnames(names1{k});
    h1(k)=uimenu(hmainplugin,'Label',names1{k});
    modulethere2=false;
    for l=1:length(names2)
        h2(k,l)=uimenu(h1(k),'Label',names2{l});
        names3=pluginnames(names1{k},names2{l});
        modulethere3=false;
        for m=1:length(names3)     
                pluginpath=pluginnames(names1{k},names2{l},names3{m});
                pname=pluginpath{4};
                ptype=pluginpath{5};
                if any(strcmp(nomenutypes,ptype))
                    continue
                end
                h3(k,l,m)=uimenu(h2(k,l),'Label',pname,'Callback',{@makeplugin,obj,{names1{k},names2{l},names3{m}}});
                modulethere3=true;  
                modulethere2=true;
                
                pout.(names1{k}).(names2{l}).(names3{m}).module={names1{k},names2{l},names3{m},pname,ptype};
        end
        if ~modulethere3
            delete(h2(k,l));
        end
    end
    if ~modulethere2
        delete(h1(k));
    end
end
%custom menu


obj.createGlobalSetting('customMenuFile','Directories','Configuration file for custom menu. Delete path and save to not have any custom menu.',struct('Style','file','String','settings/custommenu.txt'))
gfile=obj.getGlobalSetting('customMenuFile');
if exist(gfile,'file')
p=readstruct(gfile,{},true);
    if ~isempty(p)
        makecustommenu(obj.handle,p,obj)
    end
end 
end

function makecustommenu(handle,p,obj)
    fn=fieldnames(p);
    fn=setdiff(fn,{'module','position','name'});
    for k=1:length(fn)
        pm=p.(fn{k});
%         module{k}=pm.module;
        if isfield(pm,'position')
            pos(k)=pm.position;
        else
            pos(k)=inf;
        end

        if isfield(pm,'name')
            name{k}=pm.name;
        else
            name{k}=(fn{k});
        end
    end
    [~,indsort]=sort(pos);
    
    for k=1:length(fn)
        phere=p.(fn{indsort(k)});
        if isfield(phere,'module')
            uimenu(handle,'Label',name{indsort(k)},'Callback',{@makeplugin,obj,phere.module});
        else
            hs=uimenu(handle,'Label',name{indsort(k)});
            makecustommenu(hs,phere,obj); 
        end
    end
        
%     if any(strcmp(fieldnames(p.(fn{1})),'module')) %last level
%         for k=1:length(fn)
%             uimenu(handle,'Label',name{indsort(k)},'Callback',{@makeplugin,obj,p.(fn{indsort(k)}).module});
%         end
%     else
%         
%         for k=1:length(fn)
%             hs=uimenu(handle,'Label',name{indsort(k)});
%             makecustommenu(hs,p.(fn{indsort(k)}),obj);        
%         end
%         
%     end
end

function makeplugin(a,b,obj,pluginpath)
if ~iscell(pluginpath)&&strcmp(pluginpath,'Workflow')
    module=interfaces.Workflow;
    module.processorgui=false;
    name='Workflow';
    p.Vrim=5;
else
    module=plugin(pluginpath{:});
    name=pluginpath{end};
    p.Vrim=100;
end
    module.handle=figure('MenuBar','none','Toolbar','none','Name',name);
    module.attachPar(obj.P);
    module.attachLocData(obj.locData);
    
    p.Xrim=10;
    module.setGuiAppearence(p)
    module.makeGui;
end

function info_callback(a,b)
msgbox('Superresolution microscopy analysis platform (SMAP). Jonas Ries, EMBL, Heidelberg')
end

function globalsettings_callback(a,b,obj)
gui.GlobalParameterSettings([],obj.P);
end

function exit_callback(a,b,obj)
delete(obj)
end

function renamewindow_callback(a,b,obj)
title=obj.handle.Name;
answ=inputdlg('new name for SMAP window','Rename window',1,{title});
if ~isempty(answ)
    obj.handle.Name=answ{1};
end
jDesktop = com.mathworks.mde.desk.MLDesktop.getInstance;
jDesktop.getMainFrame.setTitle(answ{1});
clear jDesktop;
end

function simplegui_callback(hmenu,b,obj)
if strcmp(hmenu.Checked,'off')
    hmenu.Checked='on';
    obj.setPar('globalGuiState','s')
else
   hmenu.Checked='off';
   obj.setPar('globalGuiState','a')
end
end