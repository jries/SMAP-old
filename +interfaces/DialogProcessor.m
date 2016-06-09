classdef DialogProcessor<interfaces.GuiModuleInterface & interfaces.LocDataInterface
    %extends GuiModuleInterface with a results window, a Process button adn
    %Info button. This is the default class for Analyzer and Processor
    %modules
    properties
        resultstabgroup;   %handle to results  
        processorgui=true; %switch. true if process button etc are to be rendered. false if called externally (for workflow)
        showresults=false; % defined state for results
        history=false;
        simplegui=false;
        guiselector=struct('position',[],'show',false);
%         moduleinfo;
    end
    properties (SetAccess = private, GetAccess = private)
       resultshandle
    end
    methods
        function obj=DialogProcessor(varargin)
            obj@interfaces.GuiModuleInterface(varargin{:})  
        end
        
        function makeGui(obj,guidef)
            %calls makeGui@GuiModuleInterface, additionally provides info
            %and process button and checkbox to show results
            if nargin==1
                guidef=obj.guidef;
            end
            if obj.processorgui
            obj.guiPar.fontsize=obj.guiPar.fontsize-1;
            obj.guiPar.FieldHeight=obj.guiPar.FieldHeight-2;
            end
     
            makeGui@interfaces.GuiModuleInterface(obj,guidef);
            
            if obj.guiselector.show
                posh=obj.handle.Position;
                pos(1:2)=posh(1:2)+posh(3:4)-[23,22];
                pos(3:4)=[20,20];
                obj.guihandles.simplegui=uicontrol(obj.handle,'Style','togglebutton','String','A','Position',pos,'Callback',@obj.simplegui_callback);
                obj.guihandles.simplegui.TooltipString='toggle between simple and advanced GUI';
            end
            
            if obj.processorgui && ~isempty(obj.handle)
                hpos=obj.handle.Position;
                vrim=obj.guiPar.Vrim;
                obj.guihandles.showresults=uicontrol(obj.handle,'Position',[obj.guiPar.FieldWidth*2, hpos(4)-vrim+20,120,20],...
                    'FontSize',obj.guiPar.fontsize,'Style','checkbox', 'String', 'Show results',...
                    'Value',obj.showresults,'Callback',{@showresults_callback,obj});
                obj.guihandles.processgo_b=uicontrol(obj.handle,'Position',[obj.guiPar.FieldWidth*3, hpos(4)-vrim+20,100,50],...
                    'Style','pushbutton','String','Run','FontSize',obj.guiPar.fontsize*1.5,'Callback',{@processgo_callback,obj});
                obj.guihandles.info=uicontrol(obj.handle,'Position',[obj.guiPar.FieldWidth*2, hpos(4)-vrim+50,100,25],...
                    'Style','pushbutton','String','Info','FontSize',obj.guiPar.fontsize,'Callback',{@info_callback,obj});
            end
        end
        function setvisibility(obj,name)
            %shows and hides the GUI. called from module selector
            %setvisibility(visible) visible='on'/'off'
            if isvalid(obj.handle)&&~isa(obj.handle.Parent,'matlab.ui.Figure')
            set(obj.handle,'Visible',name);
            end
        end
        
        function makeResultsWindow(obj)
            %creates a window with results tabs
            obj.resultshandle=figure;
            obj.resultshandle.Visible='off';
            htab=uitabgroup(obj.resultshandle);
            obj.guihandles.resultstabgroup=htab;
            obj.resultstabgroup=obj.guihandles.resultstabgroup;
        end
        function ax=initaxis(obj,varargin)
            %initializes axis in results window
            ax=initaxis(obj.resultstabgroup,varargin{:});
        end
        function processgo(obj)
            %provides external access to run module (usually via process
            %button)
            processgo_callback(0,0,obj);
        end
        function addhistory(obj)
            p.parameters=obj.getGuiParameters(true,true);
            p.name=class(obj);
            obj.locData.addhistory(p);
        end      
        function simplegui_callback(obj,a,b)
            simplegui=obj.getSingleGuiParameter('simplegui');
%             if ~simplegui
%                 obj.guihandles.simplegui.String='A';
%             else
%                 obj.guihandles.simplegui.String='S';
%             end
            obj.fieldvisibility('guistate',simplegui);
%             obj.setguistate;
        end
        function fieldvisibility(obj,varargin)
            p=fieldvisibiltyparser(varargin);
            pard=obj.guidef;
            fn=fieldnames(pard);
            oldstate=obj.simplegui;
            switch p.guistate
                case {'S','s','simple',1,true}
                    obj.simplegui=true;
                    obj.guihandles.simplegui.String='S';
                    obj.guihandles.simplegui.Value=1;
                    
                    for k=1:length(fn)
                        if isfield(pard.(fn{k}),'Optional')&&pard.(fn{k}).Optional==true
                            if oldstate==false %if was advanced
                                obj.guihandles.(fn{k}).UserData.visibleSMAP=obj.guihandles.(fn{k}).Visible;
                            end
                            obj.guihandles.(fn{k}).Visible='off';
                        end
                    end
%                     if previeously advanced: write visible status into
%                     guihandles.field.smapVisible for optional parameters,
%                     switch off those parameters
                case {'A','a','advanced',0,false}
%                     switch on all optional fields with
%                     smapVisible=visible
                    obj.simplegui=false;
                    obj.guihandles.simplegui.String='A';
                    obj.guihandles.simplegui.Value=0;
                    for k=1:length(fn)
                        if isfield(pard.(fn{k}),'Optional')&&pard.(fn{k}).Optional==true
                            if isfield(obj.guihandles.(fn{k}).UserData,'visibleSMAP')       
                                 obj.guihandles.(fn{k}).Visible=obj.guihandles.(fn{k}).UserData.visibleSMAP;
                            else
                                obj.guihandles.(fn{k}).Visible='on';
                            end
                            
                        end
                    end

            end
            if ~iscell(p.on)
                p.on={p.on};
            end
            for k=1:length(p.on)
                if isfield(obj.guihandles,p.on{k})
                    if ~obj.simplegui||~(isfield(pard.(p.on{k}),'Optional')&&pard.(p.on{k}).Optional==true) %not an optional parameter
                        obj.guihandles.(p.on{k}).Visible='on';
                    end
                    obj.guihandles.(p.on{k}).UserData.visibleSMAP='on';
                end
            end
            if ~iscell(p.off)
                p.off={p.off};
            end            
            for k=1:length(p.off)
                if isfield(obj.guihandles,p.off{k})
                    
                    obj.guihandles.(p.off{k}).Visible='off';
                    
                    obj.guihandles.(p.off{k}).UserData.visibleSMAP='off';
                end
            end            
%             for all p.on: smapvisible='on'. if optional: switch on if state=advanced. 
%             for all p.off: smapvisible, visible='off';
        end
        function setguistate(obj,state)
            if nargin>1
                obj.simplegui=state;
            end
            if ~obj.simplegui
                state='on';
            else
                state='off';
            end
            pard=obj.guidef;
            fn=fieldnames(pard);
            for k=1:length(fn)
                if isfield(pard.(fn{k}),'Optional')&&pard.(fn{k}).Optional==true
                    obj.guihandles.(fn{k}).Visible=state;
                end
            end
        end

    end
    methods (Access=private)
%         function outhandle=addresultstab(obj,name)
%             %adds a tab to results figure. 
%             outhandle=uitab(obj.guihandles.resultstabgroup,'Title',name);
%         end
    end
end

function processgo_callback(a,b,obj)
% notify(obj.locData,'undo',recgui.simpleEvent('backup'));

obj.status(['executing ' class(obj)])
drawnow;
if isempty(obj.resultshandle)||~isvalid(obj.resultshandle)
    obj.makeResultsWindow;
end


p=obj.getAllParameters;

 p.obj=obj;
p.resultstabgroup=obj.guihandles.resultstabgroup;
if obj.processorgui
obj.resultshandle.Visible=onoff(p.showresults);
end
if isempty(obj.locData.loc)
    warning('no localization data present')
end
    results=obj.run(p);
if ~isempty(results)
    obj.setAutoResults(obj.pluginpath,results);
    if isfield(results,'clipboard')
        cl=results.clipboard;
        if ~iscell(cl)
            cl={cl};
        end
        
        for k=1:length(cl)
            if isnumeric(cl{k})
                cl{k}=num2str(cl{k});
            end
        end
        ct=sprintf('%s\t',cl{:});
        clipboard('copy',ct);
        
    end
end
if obj.history
    obj.addhistory;
end
obj.resultshandle.Visible=onoff(p.showresults);
if ~isfield(results,'error')||isempty(results.error)
    obj.status([class(obj) ' finished'])
else
    obj.status(['ERROR in ' class(obj) '. ' results.error])
end
end


function showresults_callback(object,data,obj)
if isempty(obj.resultshandle)||~isvalid(obj.resultshandle)
    obj.makeResultsWindow;
end
switch object.Value
    case 1
        state='on';
    case 0 
        state='off';
end
set(obj.resultshandle,'Visible',state)
end

function info_callback(a,b,obj)
obj.guihandles.showresults.Value=1;
showresults_callback(obj.guihandles.showresults,0,obj)
ax=obj.initaxis('Info');
hp=ax.Parent;
 htxt=uicontrol(hp,'Style','text','Units','normalized','Position',[0,0,.95,1],...
     'FontSize',obj.guiPar.fontsize,'HorizontalAlignment','left');
txt=textwrap(htxt,{obj.info.description});
 htxt.String=txt;
  htxt.Position=[0 0 1 1];
end


function pres=fieldvisibiltyparser(args)
% fields{end+1}='all';
p = inputParser;   
p.KeepUnmatched=true;
addParameter(p,'guistate','nd',@(x) any(myvalidatestring(x,{'simple','advanced','s','a','S','A','nd'})));
addParameter(p,'off',{});
addParameter(p,'on',{});

parse(p,args{:});
pres=p.Results;

end