classdef volume3D<interfaces.DialogProcessor
    properties
        zold
    end
    methods
        function obj=volume3D(varargin)        
            obj@interfaces.DialogProcessor(varargin) ;
            obj.inputParameters={'mingaussnm','mingausspix','gaussfac','pixrec','znm_min','znm_max'};
             obj.zold.changed=0;

        end
        function makeGui(obj)
            makeGui@interfaces.DialogProcessor(obj);
            set(obj.guihandles.preview,'Callback',{@preview,obj})
        end
        function out=run(obj,p)
            make3Dvolume(obj,p)           
            out=0;
        end
        function pard=pardef(obj)
            pard=pardef;
        end
        

    end
    methods(Static)
        function info=info(obj)
            info.name='2. Volume3D';
            info.class=@volume3D;
            info.tag='volume3D';
%             obj.info=info;
        end

    end
end


function make3Dvolume(obj,p)
global reconstructgui_abort
%store old z values
firstlayer=find(obj.locData.parameters.sr_layerson);
pld=obj.locData.layer(firstlayer(1)).parameters;
zold=[pld.znm_min pld.znm_max];
if p.pixauto
    pixz=p.pixrec;
else
    pixz=p.pixrecset;
end

if ~isfield(obj.zold,'value')
obj.zold.value=zold;
% obj.zold.auto=0
end
if obj.zold.changed
    if p.preview
        mz=(p.zmin+p.zmax)/2;
        zrange=[mz mz+pixz];
        obj.zold.value=zold;
        notify(obj.locData,'synchronizeGui',recgui.SynchronizeGuiData(-1,'znm',num2str(zrange(1)),num2str(zrange(2)),0));
    else
        zold=obj.zold.value;
       notify(obj.locData,'synchronizeGui',recgui.SynchronizeGuiData(-1,'znm',num2str(zold(1)),num2str(zold(2)),0));
    end
    obj.zold.changed=0;
else

if p.preview  
    mz=(p.zmin+p.zmax)/2;
    zrange=[mz mz+pixz];
else
    zrange=p.zmin:pixz:p.zmax;
end
imsize=size(obj.locData.guiData.imageshow);
imstack=zeros(imsize(1),imsize(2),imsize(3),length(zrange)-1);
for k=1:length(zrange)-1
    if reconstructgui_abort
        error('stopped manually')

    end
    notify(obj.locData,'synchronizeGui',recgui.SynchronizeGuiData(-1,'znm',num2str(zrange(k)),num2str(zrange(k+1))));
    notify(obj.locData,'renderImage');
    imstack(:,:,:,k)=obj.locData.guiData.imageshow;
end
if ~p.preview
notify(obj.locData,'synchronizeGui',recgui.SynchronizeGuiData(-1,'znm',num2str(zold(1)),num2str(zold(2))));
end

if p.filterz>0
    sigma = p.filterz;
    sizef = 3*ceil(sigma);
    x = linspace(-sizef / 2, sizef / 2, sizef);
    gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
    gaussFilter = gaussFilter / sum (gaussFilter); % normalize
    imstack=filter(gaussFilter,1,imstack,[],4);
end

if ~p.pixauto
    scale=p.pixrecset/p.pixrec;
    nfnew=round(length(zrange)-1)*scale;
    imscaled=zeros(imsize(1),imsize(2),imsize(3),nfnew);
    for k=1:imsize(3);
        imscaled(:,:,k,:)=imresize3d(squeeze(imstack(:,:,k,:)),[],[imsize(1),imsize(2),nfnew],'cubic');
    end
    imstack=imscaled;
end


if p.preview
else 
    if p.d3_color
        options.color=true;
    else
        options.color=false;
        imstack=squeeze(sum(imstack,3));
    end
    options.message=true;
    options.comp='lzw';

    [p,f]=fileparts(obj.locData.files.file(1).name);

    [f,p]=uiputfile([p filesep f '_3Dvolume.tif']);
    fileout=[p f];
    if f
    imout=uint8(imstack/max(imstack(:))*(2^8-1));
    saveastiff(imout,fileout,options)
    end
end  
end
end

function preview(a,b,obj)
p=obj.getGuiParameters.par;
p=copyfields(p,obj.parameters);
obj.zold.changed=true;
obj.run(p);
end

function pard=pardef
pard.d3_color.object=struct('Style','checkbox','String','render in color');
pard.d3_color.position=[2,1];

pard.text2.object=struct('String','zmin','Style','text');
pard.text2.position=[2,1];
pard.text3.object=struct('String','zmax','Style','text');
pard.text3.position=[3,1];

pard.zmin.object=struct('Style','edit','String',-400); 
pard.zmin.position=[2,2];
pard.zmax.object=struct('Style','edit','String',400); 
pard.zmax.position=[3,2];

pard.d3_color.object=struct('Style','checkbox','String','render in color');
pard.d3_color.position=[5,1];


pard.pixauto.object=struct('Style','checkbox','String','pixelsize in z auto','Value',1);
pard.pixauto.position=[4,1];
pard.pixrecset.object=struct('Style','edit','String','5'); 
pard.pixrecset.position=[4,2];

pard.text4.object=struct('String','filter  \sigma in z (planes)','Style','text');
pard.text4.position=[6,1];
pard.filterz.object=struct('Style','edit','String','1'); 
pard.filterz.position=[6,2];


pard.preview.object=struct('Style','checkbox','String','preview');
pard.preview.position=[2,4];

end
