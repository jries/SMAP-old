classdef volume3D<interfaces.DialogProcessor
    properties

    end
    methods
        function obj=volume3D(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.inputParameters={'sr_pixrec','linewidth_roi','layer1_'};
        end
        function makeGui(obj)
            makeGui@interfaces.DialogProcessor(obj);
        end
        function out=run(obj,p)
            make3Dvolume(obj,p)           
            out=0;
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end

function make3Dvolume(obj,p)
global SMAP_stopnow
lochere=obj.locData.copy;
[locsout,indout,hroi]=lochere.getloc({'xnm','ynm','znm','xnmline','ynmline'},'position','roi');
lochere.removelocs(~indout);


if ~isempty(locsout.xnmline)
    lochere.loc.xnm=locsout.xnmline;
    lochere.loc.ynm=locsout.ynmline;
    phere.sr_pos=[0 0];
    pos=hroi.getPosition;
    dx=(pos(2,1)-pos(1,1))*1000;
    dy=(pos(2,2)-pos(1,2))*1000;
    len=sqrt(dx^2+dy^2);
   
    phere.sr_size=[len p.linewidth_roi]/2;
elseif isvalid(hroi)
    pos=hroi.getPosition;
    phere.sr_pos=[pos(1)+pos(3)/2,pos(2)+pos(4)/2]*1000;
    phere.sr_size=[pos(3) pos(4)]*1000/2;
end

if p.pixzauto
    pixz=p.pixzrecset;
else
    pixz=p.sr_pixrec;
end

if p.pixxyauto
    pixxy=p.pixxyrecset;
else
    pixxy=p.sr_pixrec;
end

lochere.regroup;
lochere.filter;

 pall=obj.getLayerParameters;
 phere.imaxtoggle=0;
 phere.imax=10000;
 phere.sr_axes=[];
 phere.sr_pixrec=pixxy;
 
 
for k=1:length(pall)
    pall{k}=copyfields(pall{k},p);
    pall{k}=copyfields(pall{k},phere);
end

imall=TotalRender(lochere,pall,{'xnm','ynm'});

sim=size(imall.image);
ax=initaxis(p.resultstabgroup,'image');

z=p.zmin:pixz:p.zmax;
imout3=zeros(sim(1),sim(2),3,length(z));

showtic=tic;
updatetime=1;
if p.sigmazauto
    sz=lochere.loc.znm*0+p.sigmaz;
else
    sz=max((lochere.loc.locprecznm)*p.layer1_.gaussfac,p.layer1_.mingausspix*pixz);%gaussfac etc
    szg=max((lochere.grouploc.locprecznm)*p.layer1_.gaussfac,p.layer1_.mingausspix*pixz);
end

for k=1:length(z)
    dz=lochere.loc.znm-z(k);
    dzg=lochere.grouploc.znm-z(k);
    lochere.loc.intensity_render=exp(-(dz.^2/2./sz.^2))/sqrt(2*pi)./sz*pixz;
  lochere.grouploc.intensity_render=exp(-(dzg.^2/2./szg.^2))/sqrt(2*pi)./szg*pixz;

    srim=TotalRender(lochere,pall,{'xnm','ynm'});
    imout3(:,:,:,k)=imout3(:,:,:,k)+srim.composite;
    imdraw=srim.composite;
    imdraw=imdraw/max(imdraw(:));
    if toc(showtic)>updatetime
        imagesc(imdraw,'Parent',ax);
        title(z(k))
         drawnow
         showtic=tic;
    end
    if SMAP_stopnow
        break
    end
end
imagesc(imall.composite/max(imall.composite(:)),'Parent',ax);
improject=permute(squeeze(sum(imout3,1)),[3 1 2]);
ax2=initaxis(p.resultstabgroup,'projection');
imagesc(improject/max(improject(:)),'Parent',ax2);

options.color=true;
options.message=true;
options.comp='lzw';

[p,f]=fileparts(obj.locData.files.file(1).name);
[f,p]=uiputfile([p filesep f '_3Dvolume.tif']);
fileout=[p f];
if f
    imout=uint8(imout3/max(imout3(:))*(2^8-1));
    saveastiff(imout,fileout,options)
end
end


% function make3Dvolume(obj,p)
% global reconstructgui_abort
% %store old z values
% firstlayer=find(obj.locData.parameters.sr_layerson);
% pld=obj.locData.layer(firstlayer(1)).parameters;
% zold=[pld.znm_min pld.znm_max];
% if p.pixauto
%     pixz=p.pixrec;
% else
%     pixz=p.pixrecset;
% end
% 
% if ~isfield(obj.zold,'value')
% obj.zold.value=zold;
% % obj.zold.auto=0
% end
% if obj.zold.changed
%     if p.preview
%         mz=(p.zmin+p.zmax)/2;
%         zrange=[mz mz+pixz];
%         obj.zold.value=zold;
%         notify(obj.locData,'synchronizeGui',recgui.SynchronizeGuiData(-1,'znm',num2str(zrange(1)),num2str(zrange(2)),0));
%     else
%         zold=obj.zold.value;
%        notify(obj.locData,'synchronizeGui',recgui.SynchronizeGuiData(-1,'znm',num2str(zold(1)),num2str(zold(2)),0));
%     end
%     obj.zold.changed=0;
% else
% 
% if p.preview  
%     mz=(p.zmin+p.zmax)/2;
%     zrange=[mz mz+pixz];
% else
%     zrange=p.zmin:pixz:p.zmax;
% end
% imsize=size(obj.locData.guiData.imageshow);
% imstack=zeros(imsize(1),imsize(2),imsize(3),length(zrange)-1);
% for k=1:length(zrange)-1
%     if reconstructgui_abort
%         error('stopped manually')
% 
%     end
%     notify(obj.locData,'synchronizeGui',recgui.SynchronizeGuiData(-1,'znm',num2str(zrange(k)),num2str(zrange(k+1))));
%     notify(obj.locData,'renderImage');
%     imstack(:,:,:,k)=obj.locData.guiData.imageshow;
% end
% if ~p.preview
% notify(obj.locData,'synchronizeGui',recgui.SynchronizeGuiData(-1,'znm',num2str(zold(1)),num2str(zold(2))));
% end
% 
% if p.filterz>0
%     sigma = p.filterz;
%     sizef = 3*ceil(sigma);
%     x = linspace(-sizef / 2, sizef / 2, sizef);
%     gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
%     gaussFilter = gaussFilter / sum (gaussFilter); % normalize
%     imstack=filter(gaussFilter,1,imstack,[],4);
% end
% 
% if ~p.pixauto
%     scale=p.pixrecset/p.pixrec;
%     nfnew=round(length(zrange)-1)*scale;
%     imscaled=zeros(imsize(1),imsize(2),imsize(3),nfnew);
%     for k=1:imsize(3);
%         imscaled(:,:,k,:)=imresize3d(squeeze(imstack(:,:,k,:)),[],[imsize(1),imsize(2),nfnew],'cubic');
%     end
%     imstack=imscaled;
% end
% 
% 
% if p.preview
% else 
%     if p.d3_color
%         options.color=true;
%     else
%         options.color=false;
%         imstack=squeeze(sum(imstack,3));
%     end
%     options.message=true;
%     options.comp='lzw';
% 
%     [p,f]=fileparts(obj.locData.files.file(1).name);
% 
%     [f,p]=uiputfile([p filesep f '_3Dvolume.tif']);
%     fileout=[p f];
%     if f
%     imout=uint8(imstack/max(imstack(:))*(2^8-1));
%     saveastiff(imout,fileout,options)
%     end
% end  
% end
% end

function preview(a,b,obj)
p=obj.getGuiParameters.par;
p=copyfields(p,obj.parameters);
obj.zold.changed=true;
obj.run(p);
end

function pard=guidef
% pard.d3_color.object=struct('Style','checkbox','String','render in color');
% pard.d3_color.position=[2,1];

pard.text2.object=struct('String','zmin','Style','text');
pard.text2.position=[2,1];
pard.text3.object=struct('String','zmax','Style','text');
pard.text3.position=[3,1];

pard.zmin.object=struct('Style','edit','String',-400); 
pard.zmin.position=[2,2.5];
pard.zmax.object=struct('Style','edit','String',400); 
pard.zmax.position=[3,2.5];
% 
% pard.d3_color.object=struct('Style','checkbox','String','render in color');
% pard.d3_color.position=[5,1];


pard.pixxyauto.object=struct('Style','checkbox','String','set pixelsize in xy (nm) to:','Value',0);
pard.pixxyauto.position=[4,1];
pard.pixxyauto.Width=1.5;
pard.pixxyrecset.object=struct('Style','edit','String','5'); 
pard.pixxyrecset.position=[4,2.5];

pard.pixzauto.object=struct('Style','checkbox','String','set pixelsize in z (nm) to:','Value',0);
pard.pixzauto.position=[5,1];
pard.pixzauto.Width=1.5;
pard.pixzrecset.object=struct('Style','edit','String','5'); 
pard.pixzrecset.position=[5,2.5];

pard.sigmazauto.object=struct('String','set sigma_z (nm) to:','Style','checkbox');
pard.sigmazauto.position=[6,1];
pard.sigmazauto.Width=1.5;
pard.sigmaz.object=struct('Style','edit','String','20'); 
pard.sigmaz.position=[6,2.5];

% 
% pard.preview.object=struct('Style','checkbox','String','preview');
% pard.preview.position=[2,4];

end
