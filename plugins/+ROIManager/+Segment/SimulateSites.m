classdef SimulateSites<interfaces.DialogProcessor&interfaces.SEProcessor
    properties
    end
    methods
        function obj=SimulateSites(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_sitefov','se_cellpixelsize','se_siteroi'};
            obj.history=true;
        end
        
        function out=run(obj,p)  
            [locst,possites]=simulatelocs(p);
%            [poslabels,possites]=getlabels(obj,p);
%            posreappear=reappear(poslabels,p.blinks);
%            p.maxframe=100000;
%            locs=locsfrompos(posreappear,p);
%            locst=singlelocs(locs);
%            figure(88);plot(locst.xnm,locst.ynm,'+',locst.xnm_gt,locst.ynm_gt,'o')
           obj.locData.addfile(['simulated_' num2str(obj.locData.files.filenumberEnd)]);
           obj.locData.addLocData(locst);
           obj.locData.sort('filenumber','frame');
           initGuiAfterLoad(obj);
           
           
           se=obj.locData.SE;
           cell=interfaces.SEsites;
           cell.pos=[mean(locst.xnm) mean(locst.ynm)];
           cell.ID=0;
           cell.info.filenumber=obj.locData.files.filenumberEnd;
           se.addCell(cell);
           for k=1:length(possites)
               thissite=interfaces.SEsites;
               thissite.pos=[possites(k).x possites(k).y];
               thissite.info.cell=cell.ID;
               thissite.info.filenumber=cell.info.filenumber;
                % thissite.cellnumber=sitepar.currentcell.number;
        %         thissite.number=sitepar.sitelist.cellnumber+1;
                se.addSite(thissite);
           end
%            try
           se.processors.preview.updateFilelist;
           se.processors.preview.updateCelllist;
           se.processors.preview.updateSitelist; 
           se.currentsite=se.sites(1);
           se.currentcell=se.cells(1);
           se.currentfile=se.files(1);
%            catch
%            end
          out=[];
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

% function locso=singlelocs(locs)
% numl=0;
% fn=fieldnames(locs(1));
% for k=1:length(locs)
%     numl=numl+length(locs(k).(fn{1}));
% end
% 
% % locso=locs(1);
% % locso.(fn{1})(1)=0;
% for k=1:length(fn)
%     locso.(fn{k})(numl,1)=locs(1).(fn{k})(1); %initialize
% end
% 
% ind=1;
% for l=1:length(locs)
%     numlh=length(locs(l).(fn{1}));
%     for k=1:length(fn)
%         locso.(fn{k})(ind:ind+numlh-1)=locs(l).(fn{k});
%     end
%     ind=ind+numlh;
% end
% end


function load_callback(a,b,obj)
f=obj.getSingleGuiParameter('coordinatefile');
[f,p]=uigetfile({'*.*';'*.tif';'*.png';'*.csv';'*.txt';'*.mat'},'Choose coordinate file',f);
if ~f
    return
end
obj.setGuiParameters(struct('coordinatefile',[p f]));
[~,~,ext]=fileparts(f);
switch ext
    case {'.txt','.csv'}
        txt='on';
        tif='off';
    case {'.tif','.png'}
        txt='off';
        tif='on';
    case '.mat'
        l=load([p f]);
        if isfield(l,'image')
            txt='off';
            tif='on';
        else
            txt='on';
            tif='off';            
        end
    case '.m'
        cf=pwd;
        cd(p)
        [~,fh]=fileparts(f);
        l=eval(fh);
        cd(cf);
        if isfield(l,'image')
            txt='off';
            tif='on';
        else
            txt='on';
            tif='off';            
        end
end

obj.guihandles.labeling_efficiency.Visible=txt;
obj.guihandles.t_labelingefficiency.Visible=txt;
obj.guihandles.tif_density.Visible=tif;
obj.guihandles.tif_numbermode.Visible=tif;
obj.guihandles.tif_imagesizet.Visible=tif;
obj.guihandles.tif_imagesize.Visible=tif;
end

% function [locs,possites]=getlabels(obj,p)
% [~,~,ext]=fileparts(p.coordinatefile);
% image=[];
% locsall=[];
% switch ext
%     case {'.txt','.csv'}
%         plocs=readtable(p.coordinatefile);
%         plocsa=table2array(plocs);
%         locsall.x=plocsa(:,1);
%         locsall.y=plocsa(:,2);
%         if size(plocsa,2)>2
%             locsall.z=plocsa(:,3);
%         end
%         locsall=copyfields(locsall,plocs,{'x','y','z'});
%     case {'.tif','.png'}
% %         locs=getlabelstiff(obj,p);
%         image=imread(p.coordinatefile);
%         img=sum(image,3)/size(image,3); %binarize
%         image=double(img)/255;
%         
%     case '.mat'
%         l=load(p.coordinatefile);
%         
%         if isfield(l,'image')
%             image=l.image;
%         else
%             locsall=copyfields([],l,{'x','y','z'});
%         end
%     case '.m'
%         cf=pwd;
%         [ph,fh]=fileparts(p.coordinatefile);
%         cd(ph)
%         l=eval(fh);
%         cd(cf);
% %         l=eval(p.coordinatefile);
%         if isfield(l,'image')
%             image=l.image;
%         else
%             locsall=copyfields([],l,{'x','y','z'});
%         end
%     otherwise
%         display('file not identified selected')
%         return
% end
% 
% distsites=p.se_sitefov;
% if length(p.numberofsites)>1
%     numberofrows=p.numberofsites(2);
%     numberofsites=p.numberofsites(1)*p.numberofsites(2);
% else
%     numberofrows=ceil(32000/p.se_sitefov);
%     numberofsites=p.numberofsites;
% end
% 
% % numeroflines=ceil(p.numberofsites/numberofrows);
% for k=numberofsites:-1:1
%     xh=mod(k-1,numberofrows);
%     yh=ceil(k/numberofrows);
%     if ~isempty(image)
%         locsh=locsfromimage(image,p);
%     else
%         locsh=labelremove(locsall,p.labeling_efficiency);
%     end
% 
%     numlocs=length(locsh.x);
%     locsh.x=reshape(locsh.x,numlocs,1);
%     locsh.y=reshape(locsh.y,numlocs,1);
%     if p.randomrot
%         angle=2*pi*rand(1);
%         [locsh.x,locsh.y]=rotcoord(locsh.x,locsh.y,angle);
%     else
%         angle=0;
%     end
%     
%     if p.randomxy 
%         dx=(rand(1)-.5)*p.randomxyd*2;dy=(rand(1)-.5)*p.randomxyd*2;
%         locsh.x=locsh.x+dx;
%         locsh.y=locsh.y+dy;
%     else
%         dx=0;
%         dy=0;
%     end
%     
%     locs(k).x=locsh.x+xh*distsites;
%     locs(k).y=locsh.y+yh*distsites;
%     locs(k).angle=angle*ones(size(locsh.x));
%     locs(k).dx_gt=dx*ones(size(locsh.x));
%     locs(k).dy_gt=dy*ones(size(locsh.x));
%         
%     possites(k).x=xh*distsites;
%     possites(k).y=yh*distsites;
% %     figure(89);plot(locs(k).x,locs(k).y,'*')
% %     waitforbuttonpress
%     
% end
% end
% 
% 
% function locs=locsfromimage(image,p)
% 
% if p.tif_numbermode.Value==1
%     density=p.tif_density;
% else
%     %calcualte density from number of locs
% %     pdensity=mean(image(:))/(p.tif_imagesize/1000)^2;
%     density=p.tif_density/mean(image(:))/(p.tif_imagesize/1000)^2;
% end
% numtot=round(density*(p.tif_imagesize/1000)^2);
% 
% x=rand(numtot,1);
% y=rand(numtot,1);
% 
% 
% xpix=ceil(x*size(image,1));
% ypix=ceil(y*size(image,2));
% linind=sub2ind(size(image),xpix,ypix);
% keep=image(linind)>rand(numtot,1);
% 
% locs.x=(x(keep)-0.5)*p.tif_imagesize;
% locs.y=(y(keep)-0.5)*p.tif_imagesize;
% % pixf=size(image,1)/(p.tif_imagesize/1000);
% % densitypixel=density/pixf^2;
% 
% end
% 
% function locs=labelremove(locin,p)
% numl=length(locin.x);
% indin=rand(numl,1)<=p;
% locs=copystructReduce(locin,indin);
% end
% 
% function locso=reappear(locs,numblinks)
% for k=length(locs):-1:1
%     locso(k)=reappeari(locs(k),numblinks);
% end
% end
% function locso=reappeari(locs,numblinks)
% fn=fieldnames(locs);
% numlocs=length(locs.(fn{1}));
% numn=round(exprnd(numblinks,numlocs,1));
% indextra=zeros(sum(numn),1);
% idx=1;
% for k=1:length(numn)
%     indextra(idx:idx+numn(k)-1)=k;
%     idx=idx+numn(k);
% end
% indout=[(1:numlocs)'; indextra];
% for k=1:length(fn)
%     locso.(fn{k})=locs.(fn{k})(indout);
% end
% end
% 
% 
% function locs=locsfrompos(locsi,p)
% for k=length(locsi):-1:1
%     locs(k)=locsfromposi(locsi(k),p);
% end
% end
% 
% function locs=locsfromposi(locsi,p)
%     numlocs=length(locsi.x);
%     phot=exprnd(p.photons,numlocs,1);
%     
%     a=100;
%     PSF=100;
%     sa=PSF+a/12;
%     phot(phot<10)=10;
%     indin=phot>=10;
%     numlocs=sum(indin);
%     %MOrtensen
%     locprecnm=sqrt(sa^2./phot.*(16/9+8*pi*sa^2*p.background^2./phot/a^2));
%     locs.phot=single(phot(indin));
%     locs.bg=single(locprecnm(indin)*0+p.background);
%     locs.locprecnm=single(locprecnm(indin));
%     locs.frame=double(ceil(rand(numlocs,1)*p.maxframe));
%     locs.xnm=single(locsi.x(indin)+randn(numlocs,1).*locprecnm(indin));
%     locs.ynm=single(locsi.y(indin)+randn(numlocs,1).*locprecnm(indin));
%     locs.xnm_gt=single(locsi.x(indin));
%     locs.ynm_gt=single(locsi.y(indin));
%     locs.dxnm_gt=single(locsi.dx_gt(indin));
%     locs.dynm_gt=single(locsi.dy_gt(indin));
%     locs.angle=single(locsi.angle(indin));
%       
% end
function pard=guidef(obj)

% 
pard.coordinatefile.object=struct('String','plugins/+ROIManager/+Segment/hidden/MakeNPCCoordinates.m','Style','edit');
pard.coordinatefile.position=[1,1];
pard.coordinatefile.Width=3;

pard.load_button.object=struct('String','Load','Style','pushbutton','Callback',{{@load_callback,obj}});
pard.load_button.position=[1,4];

pard.tif_numbermode.object=struct('String',{{'Density (labels/um^2)','Number of labels'}},'Style','popupmenu');
pard.tif_numbermode.Width=1.5;
pard.tif_numbermode.position=[2,1];

pard.tif_density.object=struct('String','100','Style','edit');
pard.tif_density.position=[2,2.5];
pard.tif_density.Width=0.5;

pard.tif_imagesizet.object=struct('String','Image width (nm)','Style','text');
pard.tif_imagesizet.Width=1.5;
pard.tif_imagesizet.position=[2,3];

pard.tif_imagesize.object=struct('String','300','Style','edit');
pard.tif_imagesize.position=[2,4.5];
pard.tif_imagesize.Width=0.5;

pard.t_labelingefficiency.object=struct('String','Labeling efficiency','Style','text');
pard.t_labelingefficiency.position=[3,1];
pard.t_labelingefficiency.Width=2;

pard.labeling_efficiency.object=struct('String','.5','Style','edit');
pard.labeling_efficiency.Width=1;
pard.labeling_efficiency.position=[3,3];

pard.t2.object=struct('String','mean re-activations','Style','text');
pard.t2.position=[4,1];
pard.t2.Width=1.5;

pard.blinks.object=struct('String','1','Style','edit');
pard.blinks.Width=.5;
pard.blinks.position=[4,2.5];


pard.t3.object=struct('String','lifetime (frames)','Style','text');
pard.t3.position=[4,3];
pard.t3.Width=1.5;

pard.lifetime.object=struct('String','1','Style','edit');
pard.lifetime.Width=.5;
pard.lifetime.position=[4,4.5];

pard.t4.object=struct('String','mean number photons','Style','text');
pard.t4.position=[5,1];
pard.t4.Width=1.5;

pard.photons.object=struct('String','1000','Style','edit');
pard.photons.Width=.5;
pard.photons.position=[5,2.5];

pard.t5.object=struct('String','BG per pixel (photons)','Style','text');
pard.t5.position=[5,3];
pard.t5.Width=1.5;

pard.background.object=struct('String','20','Style','edit');
pard.background.Width=.5;
pard.background.position=[5,4.5];



pard.t6.object=struct('String','Number of sites','Style','text');
pard.t6.position=[8,1];
pard.t6.Width=1.5;

pard.randomrot.object=struct('String','Random rotation','Style','checkbox');
pard.randomrot.Width=1;
pard.randomrot.position=[7,1];

pard.randomxy.object=struct('String','Random position (nm):','Style','checkbox');
pard.randomxy.Width=1;
pard.randomxy.position=[7,2];

pard.randomxyd.object=struct('String','20','Style','edit');
pard.randomxyd.Width=1;
pard.randomxyd.position=[7,3];

pard.numberofsites.object=struct('String','5 5','Style','edit');
pard.numberofsites.Width=.5;
pard.numberofsites.position=[8,2.5];


pard.t7.object=struct('String','number of frames','Style','text');
pard.t7.position=[8,3];
pard.t7.Width=1.5;

pard.maxframes.object=struct('String','10000','Style','edit');
pard.maxframes.Width=.5;
pard.maxframes.position=[8,4.5];


pard.plugininfo.type='ROI_Analyze';
end
