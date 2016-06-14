classdef GuiModuleInterface<interfaces.GuiParameterInterface
    %interface to provide GUI functionality to modules. 
    %The GUI is described with a simple syntax, its appearence is freely
    %configurable
    % extends GuiParameterInterface to deal with GUI specific parameters
    % which are connected to uicontrols
    %distinction between GuiParameterInterface and GuiModuleInterface not
    %very clear, could have been one class.
    properties
        handle %figure or panel where the GUI is created
        guiPar %structure which defines global parameters for gui appearence (e.g. font size)
        
        %modules created inside a module should be put here. 
        %Many functions recursively search children to e.g. retrieve and
        %restore all gui settings or change the appearence
        children
        plugininfo %space to put info about plugin. 
        excludeFromSave={};% uicontrols with this name are not saved or read out by getGuiParameters
        propertiesToSave={}; %{fileds} obj.(field) is saved/restored with getGuiParameters.
        pluginpath %{path1,path2,filename}: where plugin was stored before creation
        guihandles %handles to gui objects are stored here. ONLY flat structure. All handles should be added, otherwise they will not be resized etc.
        simplegui=false;
        guiselector=struct('position',[],'show',false);
    end
    methods
         function obj=GuiModuleInterface(varargin)
             %varargin= obj.handle, obj.P (global parameter object)
            if nargin>0
                obj.handle=varargin{1};
            end
            if nargin>1
                obj.attachPar(varargin{2});
            else
                obj.attachPar(interfaces.ParameterData);
            end
            
            %PC-Mac differences
            if ispc
                 guiPar.fontsize=11;
                 guiPar.FieldHeight=26;
                 guiPar.tabsize1=[0    -1  546 342];
                 guiPar.tabsize2=[0    -1  540 311];    
                 guiPar.Vsep=3;
                 guiPar.Xrim=3;
                 guiPar.Vrim=2;  
            else
                guiPar.fontsize=15;
                guiPar.FieldHeight=25;
                guiPar.tabsize1=[1    0  527 323];
                guiPar.tabsize2=[0    -1  521 291];
                guiPar.Vsep=1;
                guiPar.Xrim=2;
                guiPar.Vrim=0; 
            end
            
            %initialze guiPar
              
            guiPar.Xsep=1;
            guiPar.Vsep=1;
            
            guiPar.Xpos=1;
            guiPar.Vpos=1;

            if ishandle(obj.handle)
                hpos=obj.handle.Position;
                guiPar.FieldWidth=(hpos(3)-3*guiPar.Xsep-2*guiPar.Xrim)/4;
            end
            obj.guiPar=guiPar;
         end
        
        function pout=getGuiParameters(obj,getchildren,onlyedit)
            %gets all GUI parameters (editable parts of uicontrols etc) in
            %a parsed format. 
            %p=getGuiParameters(getchildren)
            %if getchildren=true (optional): get paraemters also from children
            %if onlyedit: only editable fields are returned.
            if nargin<2
                getchildren=false;
            end
            if nargin<3
                onlyedit=false;
            end
            h=obj.guihandles;
            p=[];
            if ~isempty(h)
                fn=fieldnames(h);
                for k=1:length(fn)
                    vh=obj.getSingleGuiParameter(fn{k},onlyedit);  
                    if ~isempty(vh)
                            p.(fn{k})=vh;             
                    end
                end
            end
            psave=obj.propertiesToSave;
            for k=1:length(psave)
                p.(psave{k})=obj.(psave{k});
            end
            p.classname=class(obj);
            if isprop(obj,'pluginpath')
            p.pluginpath=obj.pluginpath;
            end
            pout=p;

            %children
            if getchildren
                if isstruct(obj.children)
                    guichildren=fieldnames(obj.children);
                    for k=1:length(guichildren)
                        ph=obj.children.(guichildren{k}).getGuiParameters(true,onlyedit);
                        if ~isempty(ph)
                            pout.children.(guichildren{k})=ph;
                        end
                    end
                end
            end        
        end
        
        function par=getSingleGuiParameter(obj,field,onlyedit)
            %get parsed value of uicontrol
            % p=getSingleGuiParameter(guifield). 
            % uicontrol needs to be stored as obj.guihandles.(field)
            if nargin <3
                onlyedit=false;
            end
            hfn=obj.guihandles.(field);
%             if onlyedit
%                 switch hfn.Style
%                     case {'edit','popupmenu','listbox','checkbox','togglebutton'}
%                         par=obj.handle2value(hfn);
%                     otherwise
%                         par=[];
%                 end
%             else
            par=obj.handle2value(hfn,onlyedit);
%             end
        end
        
        
        function fieldvisibility(obj,varargin)
            %sets the gui state (simple/advanced) and sets visibility of
            %certain fields
            %Arguments:
            %'guistate' 'simple'/'advanced' 1/0
            %'on', {fields}
            %'off', {fields}
            if nargin<2
                p.on={};p.off={};p.guistate=obj.simplegui;
            else
            p=fieldvisibiltyparser(varargin);
            end
            pard=obj.guidef;
            if ~isstruct(pard)
                return
            end
            fn=fieldnames(pard);
            oldstate=obj.simplegui;
            switch p.guistate
                case {'S','s','simple',1,true}
                    obj.simplegui=true;
                    obj.guihandles.simplegui.String='v';
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
                    obj.guihandles.simplegui.String='-';
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
        function setGuiParameters(obj,p,setchildren)
            %sets parameters in GUI uicontrols
            %setGuiParameters(p,setchildren)
            % if setchildren=true (optional): set also gui parameters in
            % children
            if isempty(p)
                return
            end
            if nargin<3
                setchildren=false;
            end
            
            if isstruct(p)
                fn=fieldnames(p);
                phere=p;
                h=obj.guihandles;
                for k=1:length(fn)
                    if isfield(h,fn{k})&&isprop(h.(fn{k}),'Style')&&~strcmp(h.(fn{k}).Style,'text')&&~any(ismember(obj.excludeFromSave,fn))
                        hs=obj.value2handle(phere.(fn{k}),h.(fn{k}));
                        h.(fn{k})=copyfields(h.(fn{k}),hs);
                    end
                end
                
                psave=obj.propertiesToSave;
                for k=1:length(psave)
                    if isfield(p,psave{k})
                        obj.(psave{k})=p.(psave{k});
                    end
                end
                %put output parameters to P
                if isstruct(obj.guihandles)
                fo=intersect(intersect(fn,obj.outputParameters),fieldnames(obj.guihandles));
                for k=1:length(fo)
                    obj.updateGuiParameter(0,0,fo{k})
                end
                end
                
            elseif iscell(p) %handle, value
                hs=obj.value2handle(p{2},p{1});
                
                copyfields(p{1},hs);
            end
            
            if setchildren&&isfield(p,'children')
                fn=fieldnames(p.children);
                for k=1:length(fn)
                    if isfield(obj.children,fn{k})
                    child=obj.children.(fn{k});
                    pchild=p.children.(fn{k});
                    child.setGuiParameters(pchild,true);
                    end
                end
            end
%             obj.initializeGuiParameters;
        end

        function p=getAllParameters(obj,inputParameters,editonly)
            % gets Gui Parameters (without children) and inputParameters
            % p=getAllParameters(inputParameters)
            % inputParameters can be omitted: then obj.inputParameters are
            % used
            if nargin<2||isempty(inputParameters)
                inputParameters=obj.inputParameters;
            end
            if nargin<3
                editonly=false;
            end
            if any(strcmp(inputParameters,'layers'))
                for k=1:obj.getPar('numberOfLayers')
                    inputParameters{end+1}=['layer' num2str(k) '_'];
                end
            end
            p=getAllParameters@interfaces.ParameterInterface(obj,inputParameters);
            p=copyfields(p,obj.getGuiParameters(false,editonly));
            
        end
        
        function p=getLayerParameters(obj,layeri,inputParameters)
            %get parameters of a specific layer as specified in input
            %Parameters
            %p=getLayerParameters(layer,inputParameters)
            % if inputParameters empty: use obj.inputParameters
            %layer is layer number or vector of layers. If layer is empty:
            %use all layers.
            %p is cell array of structures
            if nargin<3
                inputParameters=fieldnames(obj.P.par);
            end
            if nargin<2
                layeri=[];
            end
            if isempty(layeri)
                layer=1:obj.getPar('numberOfLayers');
            else
                layer=layeri;
            end
            pall=obj.getAllParameters(inputParameters,false);
            for k=1:length(layer)
                p{k}=pall;
                lp=obj.getPar('','layer',layer(k));
                p{k}=copyfields(p{k}, lp);
%                 p{k}=copyfields(p{k}, lp.rec_addpar);
            end
            if length(layeri)==1
                p=p{1};
            end
        end

        function resize(obj,factor)
            % resizes all GUI and children GUIs
            % usually called from figure.SizeChangeCallback (or similar)
            
            if ~isempty(obj.guihandles)
            fn=fieldnames(obj.guihandles);
            for k=1:length(fn)
                try
                obj.guihandles.(fn{k}).FontSize=obj.guihandles.(fn{k}).FontSize*factor;
                catch err
                end
                try
                    if strcmpi(obj.guihandles.(fn{k}).Units,'pixels')
                        obj.guihandles.(fn{k}).Position=obj.guihandles.(fn{k}).Position*factor;
                    end
                catch err
                end
                try
                    if isa(obj.guihandles.(fn{k}),'matlab.ui.control.Table')
                    obj.guihandles.(fn{k}).ColumnWidth=num2cell([obj.guihandles.(fn{k}).ColumnWidth{:}]*factor);
                    end
                catch err
                end
            end
            end
            if ~isempty(obj.children)
                ch=fieldnames(obj.children);
                for k=1:length(ch)
                    obj.children.(ch{k}).resize(factor);
                end
            end    
        end
        
        function  adjusttabgroup(obj,htg)
            %adjusts width of second tabgroup on mac
            if ispc
            else
                htg.Units='pixel';
                htg.Position(1)=htg.Position(1)-8;
                htg.Position(3)=htg.Position(3)+16;
                htg.Position(2)=htg.Position(2)-12;
                htg.Position(4)=htg.Position(4)+16;
                htg.Units='normalized';
            end
        end
        
        function setGuiAppearence(obj,p)
            % sets guiPar: (e.g. font size, field height etc)
            %setGuiAppearence(p): p structure with any of
            % Vrim=0; vertical rim (space above and below)
            % Xrim=0; horizontal rim (space right and left)
            % Xsep=1; horizontal space between controls
            % Vsep=1; vertical space between controls
            
            % Xpos=1; horizontal position of GUI (in GUI-units)
            % Vpos=1; vertical GUI position
            % fontsize 
            % FieldWidth: x-extension of 1 GUI unit (minus Xrim) in pixels
            %usually not set but calculated from size of
            % figure/panel in obj.handle
            % FieldHeight: y-extenson of 1 GUI unit: height of uicontrol, 
            
            fn=fieldnames(p);
            for k=1:length(fn)
                obj.guiPar.(fn{k})=p.(fn{k});
            end        
        end 

       function makeGui(obj,guidef)
           % renders the GUI according to guidef, then calls obj.initGui.
           % if guidef not passed on: calls obj.guidef (that is the usual
           % way of defining a GUI)
            if nargin==1
                guidef=obj.guidef;
            end
            
            if ~isempty(obj.handle)
                obj.handle.Units='pixels';
                hpos=obj.handle.Position;
                obj.guiPar.FieldWidth=(hpos(3)-3*obj.guiPar.Xsep-2*obj.guiPar.Xrim)/4;
            end
            if isstruct(guidef)
                allFields=fieldnames(guidef);
                guiPar=obj.guiPar;

                anyoptional=false;
                for k=1:length(allFields) 
                    thisField=guidef.(allFields{k});
                    if strcmp(allFields{k},'syncParameters')
                        obj.syncParameters=thisField;
                    elseif strcmp(allFields{k},'inputParameters')
                        obj.inputParameters=thisField;
                    elseif strcmp(allFields{k},'outputParameters')
                        obj.outputParameters=thisField;
                    elseif strcmp(allFields{k},'plugininfo')
                        obj.plugininfo=thisField;
                    elseif isstruct(thisField) && ~isempty(obj.handle) %results name
                        if ~isfield(thisField,'object') || ~isfield(thisField.object,'Style')
                            allFields{k}
                            thisField
                            str=['guidef definition is incomplete in classe:' class(obj)];
                            error(str)
                            
                        end
                        if isfield(thisField,'Optional')&&thisField.Optional
                            anyoptional=true;
                        end
                        h=thisField.object;
%                         h=uicontrol(obj.handle,thisField.object);
                        h.FontSize=guiPar.fontsize;
                        h.Units='pixels';
%                         set(h,'Units','pixels')

                        if isfield(thisField,'Width')
                            widthf=thisField.Width;
                        else
                            widthf=1;
                        end

                        if isfield(thisField,'Height')
                            heightf=thisField.Height;
                        else
                            heightf=1;
                        end

                        switch h.Style
                            case {'pushbutton','togglebutton'}
                                hadjust=4;
                            case 'popupmenu'
                                hadjust=0;
                            otherwise
                                hadjust=0;
                        end

                        h.Position=[(guiPar.FieldWidth)*(thisField.position(2)-1+guiPar.Xpos-1)+guiPar.Xrim, ...
                            hpos(4)-(thisField.position(1)+guiPar.Vpos-1)*guiPar.FieldHeight-guiPar.Vrim-hadjust/2-guiPar.Vsep-3, ...
                           guiPar.FieldWidth*widthf-2*guiPar.Xsep,...
                            guiPar.FieldHeight*heightf-guiPar.Vsep+hadjust];
                        
                        
                        if strcmpi(h.Style,'text')
                            h.HorizontalAlignment='left';
                        end
                        
                        hg=uicontrol(obj.handle,h);
                        
                        obj.guihandles.(allFields{k})=hg;
                        thisField=myrmfield(thisField,{'Width','Height','position','object','load'});
                        remaining=fieldnames(thisField);
%                         remaining=setdiff(fieldnames(thisField),{'Width','Height','position','object','load'});
                        for kr=1:length(remaining)
                            if isprop(hg,remaining{kr})
                                hg.(remaining{kr})=thisField.(remaining{kr});
                            end
                        end
%                         if isfield(thisField,'TooltipString')
%                             h.TooltipString=thisField.TooltipString;
%                         end
                    end
                    if anyoptional
                        obj.addSynchronization('globalGuiState',[],[],@obj.setglobalguistate);
                    end
                       
                end 
            
            end
           
            obj.initGui;
            obj.setSyncParameters;
            obj.initializeGuiParameters;
       end
        
       function initGui(obj) 
           % implement if needed in subclass. Called after makeGui
           %dummy, in case not implemented
       end
       
       function status(obj,txt)
           %set the status in the status line
           numchar=65;
           if length(txt)>numchar
               txt2{1}=txt(1:numchar);
               txt2{2}=txt(numchar+1:end);
           else
               txt2=txt;
           end
           obj.setPar('status',txt2);
       end
       function setglobalguistate(obj,a,b)
           simplestate=obj.getPar('globalGuiState');
               obj.fieldvisibility('guistate',simplestate)
       end
       function p=guidef(obj)
           %overwrite in module when defining a GUI.
           %p.fieldname.object=struct('Style','edit','String','x',...)
           %    defines uicontrol with parameters from struct
           %p.fieldname.position=[row,column] defines position of uicontrol
           %    in GUI coordinates (usually four units horizontally, and one
           %    unit vertically given by obj.guiPar.FieldHeight. These need not
           %    be integer
           %p.fieldname.Width, p.fieldname.Height: Widht and Height in GUI
           %    coordinates
           %p.fieldname.TooltipString : add tooltip
           %p.syncParameters={'parameterName',field,syncmode}, field is
           %    position of uicontrol as in obj.guihandles.(field). Adds
           %    synchronization of uicontrol with global parameters via
           %    interfaces.GuiParameterInterface.addSynchronization
           %p.inputParamters: can be defined here as cell array of char
           %p.outputParameters: same
           %p.plugininfo.name, .description: longer text which describes the module
           
           p=[];
       end
       function info=info(obj)
           %returns an info structure. With defaults if empty.
            info=obj.plugininfo;
            if isempty(info)||~isfield(info,'name')
               
                name=obj.pluginpath{end};
                [~,file]=fileparts(name);
                if ~isempty(file)
                    name=file;
                end
                info.name=name;
            end
            if ~isfield(info,'description')
                try
                [infopath,infofile]=fileparts(plugincell2path(obj.pluginpath));
                descfile=[infopath filesep infofile '.info'];
%                 if exist(descfile,'file')
                    fid=fopen(descfile);
                    idx=1;
                    line=fgetl(fid);
                    txt{idx}=line;
                    while ischar(line)
                        line=fgetl(fid);
                        idx=idx+1;
                        txt{idx}=line;
                    end
                    fclose(fid);
                    info.description=txt;
%                 else
%                     info.description=info.name;
%                 end
                catch
                    info.description=info.name;
                end
            end
            obj.plugininfo=info;
       end

    end
    
    methods (Access=private)
        function setSyncParameters(obj)
               for k=1:length(obj.syncParameters)
                   sp=obj.syncParameters{k};
                   if ~isempty(sp{2})&&ischar(sp{2})
                   h=obj.guihandles.(sp{2});
                   obj.addSynchronization(sp{1},h,sp{3:end});
                   end
               end

               for k=1:length(obj.outputParameters)
                   po=obj.outputParameters{k};
%                    if ~isfield(obj.P.par,po) %not already attached
                       if isfield(obj.guihandles,po)&& isempty(obj.guihandles.(po).Callback) %is part of gui and has no callback
                           obj.guihandles.(po).Callback={@obj.updateGuiParameter,po};
                           obj.updateGuiParameter(0,0,po);
                       end     
%                    end
               end

        end
    
       function initializeGuiParameters(obj)
           po=obj.outputParameters;
           for k=1:length(po)
               if isfield(obj.guihandles,po{k})
               obj.updateGuiParameter(0,0,po{k})
               end
           end 
           pi=obj.inputParameters;
           for k=1:length(pi)
               if ~isempty(obj.guihandles)&&isfield(obj.guihandles,pi{k})
                   v=obj.getPar(pi{k});
                   obj.setGuiParameters(struct(pi{k},v))
               end
           end 
           
           ps=obj.syncParameters;
           for k=1:length(ps)
               hn=ps{k}{2};
               if ishandle(hn)
                   v=obj.getPar(ps{k}{1});
                   obj.setGuiParameters({hn,v})                   
               elseif isfield(obj.guihandles,hn)
                   v=obj.getPar(ps{k}{1});
                   obj.setGuiParameters(struct(ps{k}{2},v))
               end
           end 
       end
       
      function updateGuiParameter(obj,a,b,field)
           %writes a single parameter from the GUI as specified by field
           %updateGuiParameter(a,b,field).
           %into the global parameter structure. Usually in a callback of
           %obj.guihandles.(field). Only needed when there is an explicit
           %callback. If no callback is defined but field is among
           %obj.outputParameters, this callback is created in
           %setsyncParameters
           %a,b are space holders, because updateGuiParameters is used as
           %callback.
           p=obj.getSingleGuiParameter(field);
           if isa(obj,'interfaces.LayerInterface')
              obj.setPar(field,p,'layer',obj.layer)
           else
           obj.setPar(field,p)
           end
       end
   end
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