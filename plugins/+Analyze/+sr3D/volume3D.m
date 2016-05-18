classdef volume3D<interfaces.DialogProcessor
    % volume3D renders dataset as 3D volumes with 3D elliptical Gaussiancs
    % corresponding to locprecnm and locprecznm.
    properties
        imagestack
    end
    methods
        function obj=volume3D(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.inputParameters={'sr_pixrec','linewidth_roi','layer1_'};
            obj.showresults=true;
        end
        function makeGui(obj)
            makeGui@interfaces.DialogProcessor(obj);
        end
        function out=run(obj,p)
            make3Dvolume(obj,p)           
            out=[];
        end
        function pard=guidef(obj)
            pard=guidef(obj);
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

 phere.sr_axes=[];
 phere.sr_pixrec=pixxy;
 
 
for k=1:length(pall)
    pall{k}=copyfields(pall{k},p);
    pall{k}=copyfields(pall{k},phere);
end
 phere.imaxtoggle=1;
  phere.imax=1;
imall=TotalRender(lochere,pall,{'xnm','ynm'});
 phere.imaxtoggle=0;
 phere.imax=imall.imax;
 
 for k=1:length(pall)
    pall{k}=copyfields(pall{k},p);
    pall{k}=copyfields(pall{k},phere);
end
sim=size(imall.image);
ax=obj.initaxis('image');

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
        axis(ax,'equal')
         drawnow
         showtic=tic;
    end
    if SMAP_stopnow
        break
    end
end
imagesc(imall.composite/max(imall.composite(:)),'Parent',ax);
improject=permute(squeeze(sum(imout3,1)),[3 1 2]);

ax2=obj.initaxis('projection');
imagesc(improject/max(improject(:)),'Parent',ax2);
axis(ax2,'equal')
obj.imagestack=imout3;
end

function save_callback(a,b,obj)
options.color=true;
options.message=true;
options.comp='lzw';
[p,f]=fileparts(obj.locData.files.file(1).name);
[f,p]=uiputfile([p filesep f '_3Dvolume.tif']);
fileout=[p f];
if f
    imout3=obj.imagestack;
    imout=uint8(imout3/max(imout3(:))*(2^8-1));
    saveastiff(imout,fileout,options)
end
end

function openfiji_callback(a,b,obj)
disp('not yet implemented. Please save and open manually in fiji')
%     imout3=obj.imagestack;
%     imout=uint8(imout3/max(imout3(:))*(2^8-1));
%     mij=openfiji(obj);
%     title=obj.getPar('layer1_').ch_filelist.selection;
%     mij.createColor(title,imout,true);
end

function pard=guidef(obj)
pard.text2.object=struct('String','zmin','Style','text');
pard.text2.position=[2,1];
pard.text3.object=struct('String','zmax','Style','text');
pard.text3.position=[3,1];

pard.zmin.object=struct('Style','edit','String','-400'); 
pard.zmin.position=[2,2.5];
pard.zmax.object=struct('Style','edit','String','400'); 
pard.zmax.position=[3,2.5];

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

pard.save.object=struct('Style','pushbutton','String','Save','Callback',{{@save_callback,obj}});
pard.save.position=[2,4];

pard.openfiji.object=struct('Style','pushbutton','String','open in Fiji','Callback',{{@openfiji_callback,obj}});
pard.openfiji.position=[3,4];


pard.plugininfo.name='3D volume';
pard.plugininfo.type='ProcessorPlugin';

pard.plugininfo.description= 'volume3D renders dataset as 3D volumes with 3D elliptical Gaussiancs corresponding to locprecnm and locprecznm.';
end
