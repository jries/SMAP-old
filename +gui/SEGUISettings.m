classdef SEGUISettings< interfaces.SEProcessor
    properties
        SEpreview
    end
    methods
        function obj=SEGUISettings(varargin)
            obj@interfaces.SEProcessor(varargin{:})    
        end
        function pard=pardef(obj)
            pard=pardef;
        end
        function initGui(obj)
            initGui@interfaces.SEProcessor(obj);
            set(obj.guihandles.showSE,'Callback',{@make_siteexplorer,obj})
            set(obj.guihandles.redrawall,'Callback',{@redrawall_callback,obj})
            set(obj.guihandles.clearall,'Callback',{@clearall_callback,obj})
            
%             addlistener(obj.SE.locData,'loaded',@obj.loaded_notify);
        end

%         function loaded_notify(obj,lb,eventdata)
%             obj.updateParameters;
%         end
    end
end

function make_siteexplorer(data,b,obj)
if data.Value
    if isempty(obj.SEpreview)||~isvalid(obj.SEpreview.handle)
        pfig=figure(205);
        delete(pfig.Children)
        SEpreview=gui.SEExploreGui(pfig,obj.P);
        obj.SE.processors.preview=SEpreview;
        
        SEpreview.attachLocData(obj.SE.locData);
        SEpreview.attachSE(obj.SE);
        SEpreview.makeGui;
        obj.SEpreview=SEpreview;
        obj.setPar('se_viewer',SEpreview);
    else
        set(obj.SEpreview.handle,'Visible','on')
    end
elseif ~isempty(obj.SEpreview)&&ishandle(obj.SEpreview.handle)
    set(obj.SEpreview.handle,'Visible','off')
end
        
end

function redrawall_callback(a,b,obj)
obj.SEpreview.redrawall;

end

function clearall_callback(a,b,obj) 
obj.SEpreview.clearall;

end

function pard=pardef

pard.text3.object=struct('String','FoV (nm)','Style','text');
pard.text3.position=[1,2];
pard.text4.object=struct('String','pixelsize (nm)','Style','text');
pard.text4.position=[1,3];
pard.text5.object=struct('String','ROI (nm)','Style','text');
pard.text5.position=[1,4];

pard.text1.object=struct('String','Site','Style','text');
pard.text1.position=[2,1];

pard.text2.object=struct('String','Cell','Style','text');
pard.text2.position=[3,1];

% pard.text6.object=struct('String','File','Style','text');
% pard.text6.position=[4,1];

pard.se_sitefov.object=struct('Style','edit','String',500); 
pard.se_sitefov.position=[2,2];
pard.se_cellfov.object=struct('Style','edit','String',5000); 
pard.se_cellfov.position=[3,2];
pard.se_sitepixelsize.object=struct('Style','edit','String',3); 
pard.se_sitepixelsize.position=[2,3];
pard.se_cellpixelsize.object=struct('Style','edit','String',10); 
pard.se_cellpixelsize.position=[3,3];
% pard.filepixelsize.object=struct('Style','edit','String',1); 
% pard.filepixelsize.position=[4,3];

pard.se_siteroi.object=struct('Style','edit','String',300); 
pard.se_siteroi.position=[2,4];

pard.se_imaxcheck.object=struct('Style','checkbox','String','Set Imax for sites to: [ch1 ch2 ...]','Value',0);
pard.se_imaxcheck.position=[5,1];
pard.se_imaxcheck.Width=2;

pard.se_imax.object=struct('Style','edit','String','10','Value',1);
pard.se_imax.position=[5,3];

pard.se_rotate.object=struct('Style','checkbox','String','rotate','Value',0);
pard.se_rotate.position=[6,1];
% pard.autouptdate.object=struct('Style','checkbox','String','auto update','Value',0);
% pard.autouptdate.position=[6,2];
pard.se_drawboxes.object=struct('Style','checkbox','String','draw boxes','Value',0);
pard.se_drawboxes.position=[6,2];

pard.redrawall.object=struct('Style','pushbutton','String','redraw all','Value',0);
pard.redrawall.position=[8.5,2];

pard.clearall.object=struct('Style','pushbutton','String','clear all','Value',0);
pard.clearall.position=[10,4];

pard.showSE.object=struct('Style','togglebutton','String','show ROI manager','Value',0);
pard.showSE.position=[9,1];
pard.showSE.Height=2;

pard.outputParameters={'se_sitefov','se_cellfov','se_sitepixelsize','se_cellpixelsize','se_siteroi','se_drawboxes','se_rotate','se_imax','se_imaxcheck'};
end