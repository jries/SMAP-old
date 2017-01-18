
classdef SEEvaluationGui< interfaces.SEProcessor
    properties
        pluginnames
        processors %later to SE
        currentmodule
        moduleselection={};
    end
    methods
        function obj=SEEvaluationGui(varargin)
            obj@interfaces.SEProcessor(varargin{:})
            if nargin>0
                obj.handle=varargin{1};  
            end        
            obj.processors{1}=interfaces.SEEvaluationProcessor;
            obj.currentmodule.number=1;
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function initGui(obj)
            initGui@interfaces.SEProcessor(obj);
            
            pos=obj.guihandles.addmodule.Position;
            pos2=obj.guihandles.removemodule.Position;
            pos(2)=pos(2)-obj.guiPar.FieldHeight*5;
            pos(3)=pos2(1)+pos2(3)-pos(1);
            pos(4)=pos(4)*4;
            
            obj.guihandles.modules=uitable(obj.handle,'Position',pos,'FontSize',obj.guiPar.fontsize,...
                'ColumnEditable',[true false],'RowName',[],'ColumnName',[],...
                'ColumnWidth',{20,pos(3)-20},'CellSelectionCallback',{@modules_callback,obj});
                
            obj.pluginnames=pluginnames('ROIManager','Evaluate');
%             obj.pluginnames=plugins.plugins('Siteexplorer','evaluator');
            
            set(obj.guihandles.addmodule,'Callback',{@addmodule_callback,obj})
            set(obj.guihandles.removemodule,'Callback',{@removemodule_callback,obj})
            obj.guihandles.preview.Callback={@preview_callback,obj};
            addmodule(obj,'generalStatistics');
            
        end 
        function evaluate(obj,site)
            evaluatesite(obj,site);
        end
        function redrawall(obj,a,b)
            obj.SE.processors.preview.redrawall(true);
        end
        
        function setGuiParameters(obj,p,setchildren)
            setGuiParameters@interfaces.SEProcessor(obj,p,setchildren);
            obj.guihandles.modules.Data={};
            
            if isfield(p,'children')
                modules=fieldnames(p.children);
                for k=1:length(modules);
                    mh=modules{k};
                    while mh(end)>'0'&&mh(end)<'9'
                        mh=mh(1:end-1);
                    end
                    par=p.children.(modules{k});
                    addmodule(obj,mh,par);
                end

                %set on / off
                for k=1:length(modules)
                    if isfield(p,'modules')&&~isempty(find(strcmpi(p.modules.Data(:,2),modules{k}),1))
                        val=p.modules.Data{k,1};
                        obj.guihandles.modules.Data{k,1}=val; 
                    end
                end   
            end
            
        end
    end
end

function addmodule_callback(a,b,obj)
newmodule=listdlg('ListString',obj.pluginnames);
for k=1:length(newmodule)
    modulename=obj.pluginnames{newmodule(k)};
    addmodule(obj,modulename)
end
end

function removemodule_callback(a,b,obj)
if ~isempty(obj.moduleselection)
removeind=obj.moduleselection(1);
removename=obj.processors{removeind}.name;
obj.processors(removeind)=[];
obj.children=myrmfield(obj.children,removename);
obj.guihandles.modules.Data(removeind,:)=[];
end
end

function addmodule(obj,modulename,p)
dx=2;
table=obj.guihandles.modules;
obj.guiPar.Vrim=0;
sizeparent=get( obj.handle,'Position');
d=table.Data;
pos=table.Position;

modulename2=modulename;
ind=2;
while sum(sum(strcmpi(d,modulename2)))
    modulename2=[modulename num2str(ind)];
    ind=ind+1;
end
    
s=size(d);
sn=s(1)+1;
d{sn,1}=true;
d{sn,2}=modulename2;
process=plugin('ROIManager','Evaluate',modulename);
if isa(process,'interfaces.SEEvaluationProcessor')
process.attachPar(obj.P);
panel=uipanel(obj.handle,'Units','pixels','Position',[pos(1)+pos(3)+dx,dx,sizeparent(3)-pos(1)-pos(3)-2*dx,sizeparent(4)-2*dx-20],'Visible','off');
obj.guihandles.(['processor' num2str(sn)])=panel;
obj.children.(modulename2)=process;
process.name=modulename2;
process.modulename=modulename;
process.number=sn;
process.handle=panel;
% process.guidef;
% process.addGuiParameters(obj.guiPar);
process.attachLocData(obj.locData);
% process.showresults=true;
process.makeGui;
if nargin>2 && ~isempty(p)%set parameters
    process.setGuiParameters(p);
end

obj.processors{sn}=process;
obj.guihandles.modules.Data=d;
end
end


function modules_callback(object,data,obj)
obj.moduleselection=data.Indices;
if isempty(data.Indices)
    s=size(object.Data);
    num=s(1);
else
num=data.Indices(1);
end
% module=object.Data{num,2}
showpanel(obj,num)
obj.currentmodule.number=num;
obj.currentmodule.name=object.Data{num,2};
end

function showpanel(obj,num)
s=length(obj.processors);
for k=1:s
   obj.guihandles.(['processor' num2str(k)]).Visible='off'; 
end
obj.guihandles.(['processor' num2str(num)]).Visible='on';
end

function preview_callback(a,b,obj)
% module=obj.processors{obj.currentmodule.number};
% info=module.evaluate(obj.SE.currentsite);
%for ttesting:

 evaluatesite(obj,obj.SE.currentsite,1)
end

function evaluatesite(obj,site,ploton)
p=obj.getAllParameters;
if nargin<3
    ploton=p.evaluateon;
end

if ploton
d=obj.guihandles.modules.Data;
s=size(d);
for k=1:s(1)
    if d{k,1} %evaluate
        module=obj.processors{k};
        module.display=p.display;
        module.evaluate(site);
    end
end
end
end


function pard=guidef(obj)

pard.addmodule.object=struct('Style','pushbutton','String','add module');
pard.addmodule.position=[1,1];
pard.addmodule.Width=.9;

pard.removemodule.object=struct('Style','pushbutton','String','remove');
pard.removemodule.position=[1,1.9];
pard.removemodule.Width=.6;

pard.preview.object=struct('Style','pushbutton','String','preview');
pard.preview.position=[11,1];

pard.redrawall.object=struct('Style','pushbutton','String','redraw all','Callback',@obj.redrawall);
pard.redrawall.position=[10,1.5];

pard.evaluateon.object=struct('Style','checkbox','String','evaluate on','Value',1);
pard.evaluateon.position=[7,1];
pard.evaluateon.Width=1.5;
pard.display.object=struct('Style','checkbox','String','display','Value', 1);
pard.display.position=[8,1];
pard.display.Width=1.5;

pard.se_keeptempimages.object=struct('Style','checkbox','String','keep temp imgs','Value', 0);
pard.se_keeptempimages.position=[9,1];
pard.se_keeptempimages.Width=1.5;

pard.outputParameters={'se_keeptempimages'};
end