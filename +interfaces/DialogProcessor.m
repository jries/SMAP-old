classdef DialogProcessor<interfaces.GuiModuleInterface & interfaces.LocDataInterface
    %extends GuiModuleInterface with a results window, a Process button adn
    %Info button. This is the default class for Analyzer and Processor
    %modules
    properties
        resultstabgroup;   %handle to results  
        processorgui=true; %switch. true if process button etc are to be rendered. false if called externally (for workflow)
        showresults=false; % defined state for results
        history=false;
%         moduleinfo;
    end
    properties (SetAccess = private, GetAccess = private)
       resultshandle
    end
    methods
        function obj=DialogProcessor(varargin)
            obj@interfaces.GuiModuleInterface(varargin{:})  
        end
        
        function makeGui(obj,pardef)
            %calls makeGui@GuiModuleInterface, additionally provides info
            %and process button and checkbox to show results
            if nargin==1
                pardef=obj.pardef;
            end
            if obj.processorgui
            obj.guiPar.fontsize=obj.guiPar.fontsize-1;
            obj.guiPar.FieldHeight=obj.guiPar.FieldHeight-2;
            end
     
            makeGui@interfaces.GuiModuleInterface(obj,pardef);
            
            if obj.processorgui && ~isempty(obj.handle)
                hpos=obj.handle.Position;
                vrim=obj.guiPar.Vrim;
                obj.guihandles.showresults=uicontrol(obj.handle,'Position',[obj.guiPar.FieldWidth*2, hpos(4)-vrim+20,120,20],...
                    'FontSize',obj.guiPar.fontsize,'Style','checkbox', 'String', 'showresults',...
                    'Value',obj.showresults,'Callback',{@showresults_callback,obj});
                obj.guihandles.processgo_b=uicontrol(obj.handle,'Position',[obj.guiPar.FieldWidth*3, hpos(4)-vrim+20,100,50],...
                    'Style','pushbutton','String','Process','FontSize',obj.guiPar.fontsize*1.5,'Callback',{@processgo_callback,obj});
                obj.guihandles.info=uicontrol(obj.handle,'Position',[obj.guiPar.FieldWidth*2, hpos(4)-vrim+50,100,25],...
                    'Style','pushbutton','String','Info','FontSize',obj.guiPar.fontsize,'Callback',{@info_callback,obj});
            end
%             if ~isempty(obj.info)
%                 obj.plugininfo{1}=obj.info.name;
%             else
%                 obj.plugininfo{1}=class(obj);
%             end  
%             if isfield(pardef,'plugininfo')  
%                 li=length(pardef.plugininfo);
%                 obj.plugininfo(2:li+1)=pardef.plugininfo;
%             else
%                 obj.plugininfo{2}='no info';
%                 
%             end
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
obj.status([class(obj) ' finished'])
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
ax=initaxis(obj.resultstabgroup,'Info');
hp=ax.Parent;
 htxt=uicontrol(hp,'Style','text','Units','normalized','Position',[0,0,1,1],...
     'FontSize',obj.guiPar.fontsize,'HorizontalAlignment','left');
htxt.String=obj.info.description;
end