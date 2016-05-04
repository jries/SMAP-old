classdef GuiChannel< interfaces.LayerInterface
    properties
        shift
        colorrange
        quantilestore=-4;
    end
    methods
        function obj=GuiChannel(varargin)
            obj@interfaces.LayerInterface(varargin{:});

            obj.guiPar.par=[];
            obj.guiPar.mincall=[];
            obj.guiPar.maxcall=[];
            obj.guiPar.srmodes={'normal','z','field'};
            obj.outputParameters={'ch_filelist','channels','layercheck','rendermode','render_colormode','renderfield',...
                'groupcheck','imaxtoggle','imax','lut','remout','shift','shiftxy_min','shiftxy_max','layer','colorrange',...
                'znm_min','znm_max','frame_min','frame_max'};

        end
        
        function pard=guidef(obj)
            pard=guidef(obj);
        end

        function makeGui(obj)
            makeGui@interfaces.GuiModuleInterface(obj);
            hmain=obj.getPar('filterpanel');
            pp=hmain.Position;
            
            %filtertable gui
            hfigf=uipanel('Parent',hmain,'Units','pixel','Position',[0 200 pp(3) 130]);
            filtergui=gui.GuiFilterTable(hfigf,obj.P);
            filtergui.layer=obj.layer;
            filtergui.attachLocData(obj.locData);
            filtergui.makeGui;
            obj.children.filterTableGui=filtergui;

            %histogram gui
            hfigh=uipanel('Parent',hmain,'Units','pixel','Position',[0 0 pp(3) 200]);

            histgui=gui.GuiHistSlider(hfigh,obj.P);
            histgui.layer=obj.layer;
            histgui.attachLocData(obj.locData);
            histgui.makeGui;
            
            obj.children.histgui=histgui;
            
            
            
               %menu to detach Format Gui
            f=getParentFigure(obj.handle);
            c=uicontextmenu(f); 
%             c2=uicontextmenu(f); 
            hfigf.UIContextMenu=c;
            hfigh.UIContextMenu=c;
            
            m4 = uimenu(c,'Label','detach','Callback',{@detach_callback,obj,hmain});
%              m5 = uimenu(c2,'Label','detach','Callback',{@detach_callback,obj,obj.handle});
           
        end
        
        function initGui(obj)
            h=obj.guihandles;
            fn=fieldnames(h);
            for k=1:length(fn)
                h.(fn{k}).Callback={@parchanged_callback,obj,fn{k}};
            end
            
            callobj=obj;
            h.ch_filelist.Callback={@callobj.filelist_callback,1};

            h.rendermode.Callback={@setvisibility,obj,'rendermode'};
            h.imaxtoggle.Callback={@imaxtoggle_callback,obj};
            h.render_colormode.Callback={@render_colormode_callback,obj}; 
            h.renderfield.Callback={@renderfield_callback,obj};  
            h.default_button.Callback=@obj.default_callback;
            h.externalrender.Callback={@renderpar_callback,obj};
            
            guimodules=obj.getPar('menu_plugins');
            fn=fieldnames(guimodules.Analyze.render);
            h.externalrender.String=[{'Select'},fn];
            
            cfields={ 'colorfield','String', 'locprecnm','znm','PSFxnm','locprecznm','frame','shiftxy'};
            for k=1:length(cfields)
                 h.([cfields{k} 'b']).Callback={@obj.updatefields_callback,cfields{k}};
                 h.([cfields{k} '_min']).Callback={@obj.updatefields_callback,cfields{k}};
                 h.([cfields{k} '_max']).Callback={@obj.updatefields_callback,cfields{k}};
            end
            obj.addSynchronization([obj.layerprefix 'layercheck'],h.layercheck,'Value',{@callobj.updatelayer_callback,'layercheck'})       
            h.parbutton.Callback={@renderpar_callback,obj};
            
%             cross-layer synch
             callobj=obj;
            obj.addSynchronization('filelist_short',obj.guihandles.ch_filelist,'String',{@callobj.filelist_callback,1})
           
            obj.addSynchronization([obj.layerprefix 'selectedField'],[],[],{@callobj.selectedField_callback})
            obj.addSynchronization('locFields',[],[],{@callobj.setcolorfield_callback})
            obj.addSynchronization([obj.layerprefix 'groupcheck'],h.groupcheck,'String')
            obj.addSynchronization([obj.layerprefix 'channels'],h.channels,'String')
            obj.addSynchronization([obj.layerprefix 'shiftxy_min'],h.shiftxy_min,'String',{@callobj.updatefields_callback,'shiftxy'})
            obj.addSynchronization([obj.layerprefix 'shiftxy_max'],h.shiftxy_max,'String',{@callobj.updatefields_callback,'shiftxy'})
            obj.addSynchronization(['frame_min'],[],[],{@callobj.updateframes_callback})
            obj.addSynchronization(['frame_max'],[],[],{@callobj.updateframes_callback})
            obj.guihandles=h;
%             obj.shift=[str2double(obj.guihandles.shiftxy_min.String), str2double(obj.guihandles.shiftxy_max.String)];
            
            recpar=renderpardialog([],1);
            p=obj.getAllParameters;
            layerp=copyfields(p,recpar);
            obj.setPar(obj.layerprefix,layerp);
            
            setvisibility(0,0,obj)
            obj.updateLayerField;
%             obj.filelist_callback;
            obj.updatelayer_callback;
        end
        
        function updateframes_callback(obj,a,b)
            p.frame_min=obj.getPar('frame_min');
            p.frame_max=obj.getPar('frame_max');
            obj.setGuiParameters(p);
        end
%         function externalrender_callback(obj,a,b)
%             renderpar_callback(0,0,obj)
%         end
        
        function filelist_callback(obj,x,y,z)
%             fn=obj.guihandles.ch_filelist.Value;
% %             if ~isempty(obj.getPar('selectedField')
%                 sf={'filenumber',fn,fn,true};
%                 obj.setPar('selectedField',sf,'layer',obj.layer)
%                 obj.setPar('selectedField',sframe,'layer',obj.layer)
%             end
                znm=obj.locData.getloc('znm');
                if isfield(znm,'znm')&&any(znm.znm~=0) % loical: dont update hist slider
                     obj.updatefields_callback(0,0,'znm',[],false)
                     obj.updatefields_callback(0,0,'locprecznm',[],false)
                     obj.updatefields_callback(0,0,'PSFxnm',false,false)
                else
                    obj.updatefields_callback(0,0,'PSFxnm',true,false)
                    obj.updatefields_callback(0,0,'znm',false,false)
                    obj.updatefields_callback(0,0,'locprecznm',false,false)
                end                   
                obj.updatefields_callback(0,0,'locprecnm',true,true)
                obj.updatefields_callback(0,0,'frame',false,false)

            obj.updateLayerField;
             setvisibility(0,0,obj)
        end
        
        function setcolorfield_callback(obj)
            if isempty(obj.locData.loc), return; end
            numlocs=length(obj.locData.loc.frame);
            obj.locData.loc.colorfield=(1:numlocs)'/numlocs*2-0.5;
            if ~isempty(obj.locData.grouploc)
                numlocs=length(obj.locData.grouploc.frame);
                obj.locData.grouploc.colorfield=(1:numlocs)'/numlocs;
            end
        end
         
        function updateLayerField(obj,field,value)
            if nargin<2 %update all
                layerp=obj.getPar(obj.layerprefix);
            	p=obj.getAllParameters;
                layerp=copyfields(layerp,p);
                obj.setPar(obj.layerprefix,layerp);

                fn=obj.guihandles.ch_filelist.Value;
%             if ~isempty(obj.getPar('selectedField')
                sf={'filenumber',fn,fn,true};
                obj.setPar('selectedField',sf,'layer',obj.layer)
                
            else
                
                if nargin<3
                    value=obj.getSingleGuiParameter(field);
                end
                %reads field from gui, writes it in P.par.layerX_.(field)
                p=obj.getPar(obj.layerprefix);
                p.(field)=value;
                obj.setPar(obj.layerprefix,p);
            end
        end
        
        function updatelayer_callback(obj,a,b,field)
            layeron=obj.getPar([obj.layerprefix 'layercheck']);
            layerp=obj.getPar(obj.layerprefix);
            layerp.layercheck=layeron;
            obj.setPar(obj.layerprefix,layerp);
        end
        
        function setlayer(obj,layer)
            if obj.layer==layer
                state='on';
                
            else
                state='off';
            end
            obj.children.filterTableGui.handle.Visible=state;
            obj.children.histgui.handle.Visible=state;
            
%             obj.filelist_callback;
        end
        
        function selectedField_callback(obj)
            sfield=obj.getPar('selectedField','layer',obj.layer);
            field=sfield{1};
            
            minf=num2str(sfield{2});
            maxf=num2str(sfield{3});
            p.([field '_min'])=minf;
            p.([field '_max'])=maxf;
            obj.setGuiParameters(p);
            
            obj.updateLayerField([field '_min'],sfield{2})
            obj.updateLayerField([field '_max'],sfield{3})
            
            allfields={'PSFxnm','znm','locprecznm','locprecnm','frame'};
            if any(strcmp(allfields,field))
                if ~isempty(sfield{4})
                    if sfield{4}
                        color=[1 1 1]*.94;
                    else
                        color=[1 1 1]*.7;
                    end
                    obj.guihandles.([field '_min']).BackgroundColor=color;
                    obj.guihandles.([field '_max']).BackgroundColor=color;
                end
            end

        end
        function updatefields_callback(obj,f1, data, field,filteron,updatehist)
            if nargin <6
                updatehist=true;
            end
              
            if nargin<4||isempty(field)
                field=f1;
                filteron=[];
            
            elseif nargin<5
                filteron=true;
            end
            obj.updateLayerField([field '_min'])
            obj.updateLayerField([field '_max'])
            if strcmp(field,'colorfield');
                filteron=false;
            end
            sf={field,obj.getSingleGuiParameter([field '_min']),obj.getSingleGuiParameter([field '_max']),filteron,updatehist};
            obj.setPar('selectedField',sf,'layer',obj.layer)
            obj.selectedField_callback;
        end
        
        function default_callback(obj,a,b)
            deffile=[ pwd '/settings/temp/Channel_default.mat'];
            fh=getParentFigure(obj.handle);
            modifiers = get(fh,'currentModifier');
            if ismember('shift',modifiers);
                %save
                defaultChannelParameters=obj.getGuiParameters;
                save(deffile,'defaultChannelParameters');
                disp('default parameters saved.');
                obj.status('default parameters saved.');
            else
                if exist(deffile,'file')
                    l=load(deffile);
                    p=l.defaultChannelParameters;
                    p=rmfield(p,{'ch_filelist','render_colormode','renderfield'});
                    obj.setGuiParameters(p);
                    obj.status('default parameters restored.');
                    %update
                    obj.updateLayerField;

                    
                else
                    disp('no default configuration found at: /settings/temp/Channel_default.mat')
                end
            end
        end
        
        function setGuiParameters(obj,varargin)
            
            setGuiParameters@interfaces.GuiModuleInterface(obj,varargin{:});  
           
            if length(varargin)>1&&varargin{2}==true
%             %update fields to histogram
                obj.updateLayerField;

            end
        end

    end
end

function imaxtoggle_callback(object,data,obj)
if(object.Value)
    object.String='quantile';
    obj.guihandles.imax.String=num2str(obj.quantilestore);
else
    object.String='Imax';
    obj.quantilestore=str2num(obj.guihandles.imax.String);
    try
    imax=obj.locData.layer(obj.layer).images.finalImages.imax;
    catch
        imax=1;
    end
    obj.guihandles.imax.String=num2str(imax);
end
obj.updateLayerField('imaxtoggle');
obj.updateLayerField('imax');
end



function render_colormode_callback(object,data,obj)
p=obj.getAllParameters;
switch p.render_colormode.selection
    case {'normal'}
        cmin=0;
        cmax=1;
        obj.setcolorfield_callback
    case 'z'
        if isfield(obj.locData.loc,'znm')
        obj.locData.loc.colorfield=obj.locData.loc.znm;
        obj.locData.grouploc.colorfield=obj.locData.grouploc.znm;
        cmin=-300;
        cmax=300;
        else 
            disp('no z information')
            return
        end
    case 'field'
        renderfield_callback(object, 0,obj)
        return
%         if ~isempty(obj.colorrange)&&length(obj.colorrange.mincall)>=p.renderfield.Value
%             cmin=obj.colorrange.mincall(p.renderfield.Value);
%             cmax=obj.colorrange.maxcall(p.renderfield.Value);
%         else
%             cmin=0;
%             cmax=1;
%         end
    otherwise %tiff file
            cmin=0;
            cmax=1;
end
% obj.setPar([obj.layerprefix 'colorfield_min'],num2str(cmin),'String')
% obj.setPar([obj.layerprefix 'colorfield_max'],num2str(cmax),'String')
obj.guihandles.colorfield_min.String=num2str(cmin);
obj.guihandles.colorfield_max.String=num2str(cmax);
obj.updateLayerField('render_colormode');
obj.updateLayerField('colorfield_min');
obj.updateLayerField('colorfield_max');
setvisibility(0,0,obj)
end

function renderfield_callback(object, handle,obj)
if isempty(obj.locData.loc)
    return
end
p=obj.getSingleGuiParameter('renderfield');
field=p.selection;
v=obj.locData.getloc(field,'layer',obj.layer).(field);
obj.locData.loc.colorfield=obj.locData.loc.(field);
obj.locData.grouploc.colorfield=obj.locData.grouploc.(field);
q=myquantilefast(v,[0.01,0.99]);
dx=10^floor(log10(q(2)/100));
minv=round(q(1)/dx)*dx;
maxv=round(q(2)/dx)*dx;            
obj.colorrange.mincall(p.Value)=minv;
obj.colorrange.maxcall(p.Value)=maxv;
obj.guihandles.colorfield_min.String=num2str(minv);
obj.guihandles.colorfield_max.String=num2str(maxv);
% render_colormode_callback(object,0,obj)      
obj.updateLayerField('render_colormode');
obj.updateLayerField('renderfield');
obj.updateLayerField('colorfield_min');
obj.updateLayerField('colorfield_max');
 setvisibility(0,0,obj)
end

function renderpar_callback(object,data,obj)
p=obj.getGuiParameters;

switch p.rendermode.selection
case 'Other'
    if strcmp(p.externalrender.selection,'Select')
        return
    end
    rendermodules=obj.getPar('rendermodules');
    
    if length(rendermodules)>=obj.layer && ~isempty(rendermodules{obj.layer})&&isvalid(rendermodules{obj.layer})
        delete(rendermodules{obj.layer}.handle)
        delete(rendermodules{obj.layer})
    end
    renderer=p.externalrender.selection;
    guiplugins=obj.getPar('menu_plugins');
    pluginpath=guiplugins.Analyze.render.(renderer).module;
    module=plugin(pluginpath{1:3});
    if strcmp(object.Style,'popupmenu')
        vis='off';
    else
        vis='on';
    end
    module.handle=figure('MenuBar','none','Toolbar','none','Name',pluginpath{end},'Visible',vis);
    module.attachPar(obj.P);
    module.attachLocData(obj.locData);
    ph.Vrim=100;
    ph.Xrim=10;
    module.setGuiAppearence(ph)
    module.makeGui;
    module.guihandles.showresults.Value=0;
    module.layer=obj.layer;
    rendermodules{obj.layer}=module;
    obj.setPar('rendermodules',rendermodules);
    
otherwise
    layerp=obj.getPar(obj.layerprefix);
    par=renderpardialog(layerp);
    if ~isempty(par)
       layerp=copyfields(layerp,par);
       obj.setPar(obj.layerprefix,layerp)       
       if par.copyall
           nl=obj.getPar('numberOfLayers');
           for k=1:nl
               layerp=obj.getPar('','layer',k);
               layerp.rec_addpar=par;
               obj.setPar('',layerp,'layer',k)
           end
       end
    end
end
end

function setvisibility(a,b,obj,field)
hgui=obj.guihandles;
if hgui.render_colormode.Value>length(hgui.render_colormode.String)
    set(hgui.render_colormode,'Value',1);
    disp('gui channel line 370')
end
p=obj.getAllParameters;


% z data?
fh={'PSFxnmb','PSFxnm_min','PSFxnm_max','znmb','znm_min','znm_max','locprecznmb','locprecznm_min','locprecznm_max',...
    'channels','text1','renderfield','groupcheck',...
    'locprecnmb','locprecnm_min','locprecnm_max'};
%tiff
if p.rendermode.Value==4 %obj.fileinfo(fileselect).istiff
    controlVisibility(hgui,fh,'off')
else
    controlVisibility(hgui,fh,'on')
    zh={'znmb','znm_min','znm_max','locprecznmb','locprecznm_min','locprecznm_max'};
    znoth={'PSFxnmb','PSFxnm_min','PSFxnm_max'};
    znm=obj.locData.getloc({'znm'},'position','all');
    
    if ~isempty(znm)&&~isempty(znm.znm)&&any(znm.znm~=0)
        controlVisibility(hgui,zh,'on')
        controlVisibility(hgui,znoth,'off')
    else
        controlVisibility(hgui,zh,'off')
        controlVisibility(hgui,znoth,'on')
    end
end

if strcmp(p.rendermode.selection,'Other')
    controlVisibility(hgui,'externalrender','on')
else
    controlVisibility(hgui,'externalrender','off')
end

%render_colormode 
hgui.render_colormode.Value=min(length(hgui.render_colormode.String),hgui.render_colormode.Value);
switch lower(p.render_colormode.selection)
case 'normal'
    controlVisibility(hgui,{'renderfield'},'off')
case 'z'
    controlVisibility(hgui,{'remout','cb','colorfield_min','colorfield_max'},'on')
    controlVisibility(hgui,{'renderfield'},'off')
case 'field'
     controlVisibility(hgui,{'renderfield','remout','cb','colorfield_min','colorfield_max'},'on')
end

% render or image
switch lower(p.rendermode.selection)
case {'hist','gauss','dl','other'}
    if isstruct(obj.locData.loc)
        set(hgui.renderfield,'String',fieldnames(obj.locData.loc));
        hgui.renderfield.Value=min(hgui.renderfield.Value,length(hgui.renderfield.String));
    end
    if hgui.render_colormode.Value>length(hgui.render_colormode.String)
        set(hgui.render_colormode,'Value',1);
        disp('gui channel line 370')
    end
    set(hgui.render_colormode,'String',{'normal','z','field'});
case {'tiff'}
    s={};
    sind=1;
    file=hgui.ch_filelist.Value;
    s{1}='empty';
    if ~isempty(obj.locData.files.file)&&isfield(obj.locData.files.file(file),'tif')
        for k=1:length(obj.locData.files.file(file).tif)
           if ~isempty(obj.locData.files.file(file).tif(k).info)
                [p,f,ext]=fileparts(obj.locData.files.file(file).tif(k).info.name);
                s{sind}=f;
                sind=sind+1;
           end
        end

        set(hgui.render_colormode,'String',s);
        set(hgui.render_colormode,'Value',min(hgui.render_colormode.Value,length(s)));
        set(hgui.colorfield_min,'String','0')
        set(hgui.colorfield_max,'String','1')
        updateLayerField(obj,'colorfield_min')
        updateLayerField(obj,'colorfield_max')
        obj.updateLayerField('render_colormode');
    end
end

if nargin>3
    updateLayerField(obj,field)
end
end




function paro=renderpardialog(par,default)
if nargin==0 || isempty(par)
    par.mingaussnm=1;
    par.mingausspix=.7;
    par.gaussfac=0.4;
    par.gamma=1;
    par.copyall=false;
end
if nargin<2||~default
[settings, button] = settingsdlg(...
    'Description', 'Render Parameters',... 
    'title'      , 'Par',...  
    {'min sigma Gauss (nm)';'mingaussnm'}, par.mingaussnm,...
    {'min sigma Gauss (pix)';'mingausspix'}, par.mingausspix,...
    {'factor Gauss';'gaussfac'}, par.gaussfac,...
    {'gamma image';'gamma'}, par.gamma,...
 {'copy to all cahnnels';'copyall'},false);

if strcmpi(button,'ok')
    par=addFields(par,settings);
    paro=par;
else
    paro=[];
end
else
    paro=par;
end
end


function parchanged_callback(object,event,obj,field)
obj.updateLayerField(field)
end


function pard=guidef(obj)
pard.layercheck.object=struct('Style','checkbox','String','','Value',1);
pard.layercheck.position=[1,1];
pard.layercheck.Width=0.2;
pard.layercheck.TooltipString='switch layer on and off';
            
pard.ch_filelist.object=struct('Style','popupmenu','String',{'File'});
pard.ch_filelist.position=[1,1.2];
pard.ch_filelist.Width=1;
pard.ch_filelist.TooltipString='which file (loc or image) to display';

pard.text1.object=struct('Style','text','String','Ch');
pard.text1.position=[1,2.2];
pard.text1.Width=0.2;

pard.channels.object=struct('Style','edit','String','0 1');
pard.channels.position=[1,2.4];
pard.channels.Width=0.6;
pard.channels.TooltipString='channels to display. Use a,b,c and a:c notation';

pard.rendermode.object=struct('Style','popupmenu','String',{{'hist','Gauss','DL','tiff','Other'}},'Value',2);
pard.rendermode.position=[1,3];
pard.rendermode.Width=1;  
pard.rendermode.TooltipString='how to render image. DL is diffraction limited';
 
pard.groupcheck.object=struct('Style','checkbox','String','group');
pard.groupcheck.position=[2,1];
pard.groupcheck.Width=.6;
pard.groupcheck.TooltipString='use grouped or ungrouped locs';

pard.imaxtoggle.object=struct('Style','togglebutton','String','quantile','Value',1);
pard.imaxtoggle.position=[2,1.6];
pard.imaxtoggle.Width=.6;
pard.imaxtoggle.TooltipString='toggle absolute intensity maximum (Imax) or quantile';

pard.imax.object=struct('Style','edit','String','-3.5');
pard.imax.position=[2,2.2];
pard.imax.Width=.6;
pard.imax.TooltipString='absolut intensity or quantile (0<q<1) or v for q=1-10^(v), v<0';

pard.render_colormode.object=struct('Style','popupmenu','String',{obj.guiPar.srmodes}); 
pard.render_colormode.position=[2,3];
pard.render_colormode.Width=1;
render_colormode.TooltipString=sprintf('normal: intensity coded, z: z color coded, \n param: select field which to color code');

pard.renderfield.object=struct('Style','popupmenu','String','field');            
pard.renderfield.position=[2,4];
pard.renderfield.Width=1;
pard.renderfield.TooltipString='field to color code';


pard.colorfieldb.object=struct('Style','pushbutton','String','c range');
pard.colorfieldb.position=[3,1];
pard.colorfieldb.Width=.6;
pard.colorfieldb.TooltipString='range of values to fill the lookup table';


pard.colorfield_min.object=struct('Style','edit','String','0');
pard.colorfield_min.position=[3,1.6];
pard.colorfield_min.Width=.6;
pard.colorfield_min.TooltipString=pard.colorfieldb.TooltipString;

pard.colorfield_max.object=struct('Style','edit','String','1');
pard.colorfield_max.position=[3,2.2];
pard.colorfield_max.Width=.6;
pard.colorfield_max.TooltipString=pard.colorfieldb.TooltipString;

% pard.lut.object=struct('Style','popupmenu','String',lutnames);

pard.lut.object=struct('Style','popupmenu','String',{mymakelut});
pard.lut.position=[3,3];
pard.lut.Width=1;
pard.lut.TooltipString='select the lookup table';

pard.remout.object=struct('Style','checkbox','String','remove out');
pard.remout.position=[3,4];
pard.remout.Width=1;
pard.remout.TooltipString=sprintf('if checked: remove loclizations outside lut. \n If unchecked: set them to maximum color');
            
pard.locprecnmb.object=struct('Style','pushbutton','String','locp');
pard.locprecnmb.position=[4,1];
pard.locprecnmb.Width=.6;


pard.locprecnm_min.object=struct('Style','edit','String','0','BackgroundColor',[1 1 1]*.7);
pard.locprecnm_min.position=[4,1.6];
pard.locprecnm_min.Width=.6;

pard.locprecnm_max.object=struct('Style','edit','String','30');  
pard.locprecnm_max.position=[4,2.2];
pard.locprecnm_max.Width=.6;

pard.znmb.object=struct('Style','pushbutton','String','z');
pard.znmb.position=[4,3.2];
pard.znmb.Width=.6;

pard.znm_min.object=struct('Style','edit','String','-500','BackgroundColor',[1 1 1]*.7);
pard.znm_min.position=[4,3.8];
pard.znm_min.Width=.6;

pard.znm_max.object=struct('Style','edit','String','500');  
pard.znm_max.position=[4,4.4];
pard.znm_max.Width=.6;
      
pard.PSFxnmb.object=struct('Style','pushbutton','String','PSF xy');
pard.PSFxnmb.position=[5,1];
pard.PSFxnmb.Width=.6;
            
pard.PSFxnm_min.object=struct('Style','edit','String','0','BackgroundColor',[1 1 1]*.7);
pard.PSFxnm_min.position=[5,1.6];
pard.PSFxnm_min.Width=.6;

pard.PSFxnm_max.object=struct('Style','edit','String','175');  
pard.PSFxnm_max.position=[5,2.2];
pard.PSFxnm_max.Width=.6;
            
pard.locprecznmb.object=struct('Style','pushbutton','String','locprec z');
pard.locprecznmb.position=[5,3.2];
pard.locprecznmb.Width=.6;

pard.locprecznm_min.object=struct('Style','edit','String','0','BackgroundColor',[1 1 1]*.7);
pard.locprecznm_min.position=[5,3.8];
pard.locprecznm_min.Width=.6;

pard.locprecznm_max.object=struct('Style','edit','String','45');  
pard.locprecznm_max.position=[5,4.4];
pard.locprecznm_max.Width=.6;
   
pard.frameb.object=struct('Style','pushbutton','String','frame');
pard.frameb.position=[6,1];
pard.frameb.Width=.6;

pard.frame_min.object=struct('Style','edit','String','0','BackgroundColor',[1 1 1]*.7);
pard.frame_min.position=[6,1.6];
pard.frame_min.Width=.6;

pard.frame_max.object=struct('Style','edit','String','inf');  
pard.frame_max.position=[6,2.2];
pard.frame_max.Width=.6;
            
pard.shiftxyb.object=struct('Style','pushbutton','String','shift x,y');
pard.shiftxyb.position=[6,3.2];
pard.shiftxyb.Width=.6;

pard.shiftxy_min.object=struct('Style','edit','String','0');
pard.shiftxy_min.position=[6,3.8];
pard.shiftxy_min.Width=.6;


pard.shiftxy_max.object=struct('Style','edit','String','0');
pard.shiftxy_max.position=[6,4.4];
pard.shiftxy_max.Width=.6;

pard.parbutton.object=struct('Style','pushbutton','String','par');
pard.parbutton.position=[1,4.5];
pard.parbutton.Width=.5;
pard.parbutton.TooltipString='Additional render paramters';

pard.externalrender.object=struct('Style','popupmenu','String','empty');
pard.externalrender.position=[1,3.9];
pard.externalrender.Width=.6;
pard.externalrender.TooltipString='External renderer';

pard.default_button.object=struct('Style','pushbutton','String','Default');
pard.default_button.position=[8,4];
pard.default_button.Width=1;
pard.default_button.TooltipString='Reset to default. Click with shift to save default.';
%%%put in again
% pard.layercolorz.object=struct('Style','checkbox','String','layers same c/z');
% pard.layercolorz.position=[7,3.8];
% pard.layercolorz.Width=1.2; 

end

function detach_callback(a,b,obj,handle)
f=figure('MenuBar','none','Toolbar','none');
handle.Parent=f;
handle.Position(1)=0;
handle.Position(2)=0;
f.Position(3:4)=handle.Position(3:4);
handle.Tag='detached';
% if strcmp(handle.Tag,'OV')
%     obj.ovdetached=true;
% end
end
