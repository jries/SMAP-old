classdef RoiAdder<interfaces.WorkflowModule
    properties
        mask
        preview
    end
    methods
        function obj=RoiAdder(varargin)
            obj@interfaces.WorkflowModule(varargin{:});
            obj.propertiesToSave={'mask'};
        end
        function pard=pardef(obj)
            pard=pardef;
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
            obj.guihandles.addroi.Callback={@addroi_callback,obj};
            obj.guihandles.removeroi.Callback={@removeroi_callback,obj};
            obj.guihandles.clearroi.Callback={@clearroi_callback,obj};
        end
        function prerun(obj,p)
            obj.preview=obj.getPar('loc_preview');
 
        end
        function resetmask(obj,a,b)
            obj.mask=[];
        end
        function dato=run(obj,data,p)

            img=data.data;%get;
            if ~isempty(img)&&all(size(obj.mask)==size(img))
                img(~obj.mask)=-1;
            end
            dato=data;%.copy;
            dato.data=img;%set(img);
            if obj.preview && data.frame==obj.getPar('loc_previewframe')
                figure(obj.getPar('loc_outputfig'));
                ax=gca;
                maxv=ax.CLim(2)/2;
                imagesc('CData',obj.mask*0+maxv,'AlphaData',double(~obj.mask)*0.5) 
            end  
        end
    end
end

function addroi_callback(a,b,obj)
p=obj.getGuiParameters;
obj.status('select roi and double clock on ROI when done...')
mask=getroi(obj.getPar('loc_outputfig'),p.roistyle.selection);
obj.status('ROI selected')
if isempty(obj.mask)||any(size(mask)~=size(obj.mask));
    obj.mask=mask;
else
    obj.mask=obj.mask|mask;
end
% figure(55);imagesc(obj.mask);
end


function removeroi_callback(a,b,obj)
p=obj.getAllParameters;
obj.status('select roi...')
mask=getroi(obj.getPar('loc_outputfig'),p.roistyle.selection);
obj.status('ROI selected')
if isempty(obj.mask)||any(size(mask)~=size(obj.mask));
    obj.mask=~mask;
else
    obj.mask=obj.mask&~mask;
end
end

function clearroi_callback(a,b,obj)
    obj.mask=[];
end

function mask=getroi(figh,roistyle)
figure(figh);
ax=gca;
switch roistyle
    case 'rectangle'
        hroi=imrect(ax);
    case 'ellipse'
        hroi=imellipse(ax);
end
hpos=wait(hroi);
him=findobj(ax.Children,'Type','Image');
k=1;
while isempty(him(k).CData)
    k=k+1;
end

mask=createMask(hroi,him(k));
end

function pard=pardef
pard.addroi.object=struct('Style','pushbutton','String','ROI to include');
pard.addroi.position=[1,1];
pard.addroi.TooltipString=sprintf('Add a region of interest in which to fit. \n Use Preview before to display image in which to select ROI. \n After selcting ROI: double click on it to confirm.');
pard.removeroi.object=struct('Style','pushbutton','String','ROI to exclude');
pard.removeroi.position=[2,1];
pard.removeroi.TooltipString=sprintf('Add a region of interest in which NOT to fit. \n Use Preview before to display image in which to select ROI. \n After selcting ROI: double click on it to confirm.');

pard.clearroi.object=struct('Style','pushbutton','String','Clear ROIs');
pard.clearroi.position=[3,1];
pard.clearroi.TooltipString=sprintf('Clear all ROIs');
pard.roistyle.object=struct('Style','popupmenu','String','rectangle|ellipse');
pard.roistyle.position=[4,1];
pard.roistyle.TooltipString=sprintf('Use rectangular or elliptical ROI.');
% pard.mask_store.object=struct('Style','text','String','');

end