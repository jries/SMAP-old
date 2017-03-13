classdef correct3Daberrations<interfaces.DialogProcessor
    properties
        SXY
    end
    methods
        function obj=correct3Daberrations(varargin)   
            obj@interfaces.DialogProcessor(varargin{:}) ;
             obj.showresults=true;
%              obj.guiPar.FieldHeight=obj.guiPar.FieldHeight-1;obj.guiPar.Vrim=obj.guiPar.Vrim-20;
        end
        function out=run(obj,p)
            p.EMon=obj.locData.files.file(1).info.EMon;
            p.RIM=obj.locData.history{1}.children.fitparamters.fitterGUI.children.MLE_GPU_Yiming.refractive_index_mismatch;
            p.dz=p.dz*p.RIM;
            out=[];
            %             * Make 3D fitting model
            % * Fit bead stacks (gel, cells) with this model 
            %     * (refractive index mismatch =1) or rescale later
            % * extract bead positions (grouping)
            disp('get beads')
            if p.beadsource.Value==1 %ROImanager
                beads=sites2beads(obj,p);
            else % segment new
                beads=segmentb(obj,p);
            end
            
            %framerange
            frange=[min(obj.locData.loc.frame) max(obj.locData.loc.frame)];
            
            % get f0 for beads
            for k=length(beads):-1:1
                   beads(k).loc.znm=beads(k).loc.znm;
                [beads(k).f0]=getf0Z(beads(k).loc,p);
            end
            f0all=([beads(:).f0]);
            badind=f0all<frange(1)|f0all>frange(2);
            beads(badind)=[];
            
            % * determine fglass, glass
            p.axhere=obj.initaxis('z0positions');
            f0glass=getf0glass(beads,p);
            p.f0glass=f0glass;
            %calculate relevant other coordinates
            
            for k=1:length(beads)
                beads(k).loc.zglass=(beads(k).loc.frame-f0glass(beads(k).filenumber))*p.dz;
                beads(k).loc.z0relative=-(beads(k).f0-beads(k).loc.frame)*p.dz;
                beads(k).loc.dzcorr=beads(k).loc.z0relative-beads(k).loc.znm;
                beads(k).f0glass=beads(k).f0-f0glass(beads(k).filenumber);
                beads(k).loc.z0glass=beads(k).f0glass*p.dz+0*beads(k).loc.zglass;
            end
            
            p.axhere=[];
            %get interpolation
            phere=p;
            phere.smoothing=[0.05 0.002];
           [ZcorrInterp]=getZinterp(beads,[],phere);
           p.cutoffrefine=300;
           [ZcorrInterp]=getZinterp(beads,ZcorrInterp,p);
           p.cutoffrefine=150;
            %calculate errors
           [err1,dzerr]=geterrors(beads,ZcorrInterp);
           err0=err1;
           goodind=find(true(length(beads),1));
           beads2=beads;
           while  length(beads2)>length(beads)/2
                cutoff=3*nanmean(err1);
                badind=(err1>cutoff|isnan(err1));
                if sum(badind)==0
                    break
                end
                goodind=goodind(~badind);
                beads2=beads(goodind);
                dzerr2=dzerr(goodind);
                [ZcorrInterp]=getZinterp(beads2,ZcorrInterp,p);
                %calculate errors
                [err1,dzerrh]=geterrors(beads2,ZcorrInterp);    
                dzerr(goodind)=dzerrh;
           end
           
           %ploto utput
           p.axhere=obj.initaxis('Interpolation');
           [ZcorrInterp]=getZinterp(beads2,ZcorrInterp,p);
           
           axhere=obj.initaxis('error');
           n=1:length(beads);
           plot(n,err0,n(goodind),err1,'*-');
           %correct beads for testing
           
           ax1=obj.initaxis('validation');
           f=ax1.Parent;
           ax2=axes(f,'Position',[0.5 0 1 1]);
           subplot(1,2,1,ax1);
           subplot(1,2,2,ax2);
       
            minv=inf;
            maxv=-inf;
           for k=1:length(beads)

               dZ=ZcorrInterp(beads(k).loc.zglass,beads(k).loc.znm);
               beads(k).loc.znmcorrected=beads(k).loc.znm+dZ;
               
               if any(goodind==k)
                   col='k.';
                   inr=abs(beads(k).loc.z0relative)<1000;
                   minv=min(min(minv,min(beads(k).loc.znm(inr))),min(beads(k).loc.znmcorrected(inr)));
                   maxv=max(max(maxv,max(beads(k).loc.znm(inr))),max(beads(k).loc.znmcorrected(inr)));
               else
                   col='r.';
               end
               plot(ax1,beads(k).loc.z0relative,beads(k).loc.znm,col)
               plot(ax2,beads(k).loc.z0relative,beads(k).loc.znmcorrected,col)
               hold(ax1,'on');
               hold(ax2,'on');
           end
           
           xlim(ax1,[-1000 1000])
           ylim(ax1,[minv maxv]);
           xlim(ax2,[-1000 1000])
           ylim(ax2,[minv maxv]);
            %get image stacks if needed
            
               % * zfitted(z0-zglass,ZObjective)
            %     * also for xfitted, yfitted, 
            %     * also spatially resolved
            % * z0-zglass=ztrue(zfitted, zObjective): interpolated or lookup table
            % * use for correction
            % * save with 3Dcal.mat
            % * fitter: instead of refractive index mismatch choose this correction.
            
            %save: either select existing 3Dcal:then it is appended. Or
            %save new (correction then with plugin)
           

        end
        
        function initGui(obj)
%             setvisible(0,0,obj)
%             beaddistribution_callback(0,0,obj)           
        end
 
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function [Zint]=getZinterp(beads,Zintold,p)
 %make big array for interpolation
 zglassall=[];z0relativeall=[];zfitall=[];idall=[];zplot=[];dzerrall=[];z0glassall=[];
            for k=1:length(beads)
                
                zglassall=double(vertcat(zglassall,beads(k).loc.zglass));
                z0glassall=double(vertcat(z0glassall,beads(k).loc.z0glass));
                z0relativeall=double(vertcat(z0relativeall,beads(k).loc.z0relative));
%                 z0all=double(vertcat(z0all,beads(k).f0*p.dz*ones(length(beads(k).loc.zglass),1)));

                zfitall=double(vertcat(zfitall,beads(k).loc.znm));
%                 fileall=double(vertcat(fileall,double(beads(k).filenumber)*ones(length(beads(k).loc.zglass),1)));
                idall=double(vertcat(idall,k*ones(length(beads(k).loc.zglass),1)));
                if 0%beads(k).f0glass<50
                    zplot=double(vertcat(zplot,beads(k).loc.dzcorr*0));
                else
                    zplot=double(vertcat(zplot,beads(k).loc.dzcorr));
                end
                
%                 if ~isempty(dzerr)
%                     dzerrall=double(vertcat(dzerrall,dzerr{k}));
%                 else
%                     dzerrall=double(vertcat(dzerrall,0*beads(k).loc.znm));
%                 end
            end
            
            zglassall=z0glassall;
%             zplot(z0glassall<150)=0;
            
            %for now don't follow up, later: remove single wrong
            %localizations.
            
            %determine range/outliers
            qzfit=myquantile(zfitall,[0.05,0.95]);
            inz=abs(z0relativeall)<800;
%             inz=inz&z0glassall>-50;
            inz=inz&(zfitall)<qzfit(2)&(zfitall)>qzfit(1);
%             figure(78);scatter3(zglassall,zfitall,z0all)
%             zplot=z0relativeall-zfitall;
            inz=inz&abs(zplot)<800;
            
            if ~isempty(Zintold)
                
                dz=Zintold(zglassall,zfitall)-zplot;

                inz=inz&abs(dz)<p.cutoffrefine;
                h=histcounts(idall(inz),(1:max(idall)+1))';
                minpoints=150;
                innump=h(idall)>minpoints;
                inz=inz&innump;
            end
            
            
            xrange=min(zglassall(inz)):100:max(zglassall(inz));
            
            yrange=[qzfit(1):10: qzfit(2)];
            [X,Y]=meshgrid(xrange,yrange);  
            %interpolation
            Z=RegularizeData3D(zglassall(inz),zfitall(inz),zplot(inz),xrange,yrange,'smoothness',p.smoothing);
             Zint=griddedInterpolant(X',Y',Z');
           
            
%              Z=RegularizeData3D(zglassall(inz),zfitall(inz),zplot(inz),xrange,yrange,'smoothness',[.1,.001]);
%             
%             Zint=griddedInterpolant(X',Y',Z');
            
            if ~isempty(p.axhere)
                scatter3(p.axhere,zglassall(inz),zfitall(inz),zplot(inz),[],idall(inz))
                xlabel(p.axhere,'zglass');ylabel(p.axhere,'zfit'); zlabel(p.axhere,'z0r');
                colormap(p.axhere,'lines')
                hold(p.axhere,'on')
                mesh(p.axhere,X,Y,Zint(X,Y),'FaceAlpha',0.2)
            end
end

function  f0glass=getf0glass(beads,p)
if isempty(p.axhere)
    f=figure;ax=gca;
else
    ax=p.axhere;
end
for k=1: max([beads(:).filenumber]) 
    indf=[beads(:).filenumber]==k;
      f0=[beads(indf).f0];
    dzh=50/p.dz;
    f0=f0(f0<dzh*40); %only look in the first 2 um
    range=0:dzh:1000;
    h=histogram(ax,f0,range);


    [mh]=max(h.Values);
    ind=find(h.Values>mh*.4,1,'first');
    f0h=range(ind);
    ind=find(f0>f0h-2*dzh&f0<f0h+2*dzh);
    f0glass(k)=mean(f0(ind));
    hold(ax, 'on');
    
end
if isempty(p.axhere)
    close(f)
else
    plot(ax,f0glass,ones(size(f0glass)),'k*')
end
end
function [err1,dzerr]=geterrors(beads,Zint)
%   figure(99)
%             hold off
            xrange=Zint.GridVectors{1};
            yrange=Zint.GridVectors{2};
  for k=1:length(beads)
        zh=double(beads(k).loc.znm);
        zglass=beads(k).loc.zglass;
        z0f=beads(k).loc.dzcorr;
        inz=abs(zh<300) & abs(z0f)<300 & (zh)<yrange(end) & (zh)>yrange(1);
%                 inz= (zh)<qzfit(2) & (zh)>qzfit(1);
        dz=Zint(zglass(inz),zh(inz))-z0f(inz);
% 
%         plot(zh(inz),dz)
%         hold on
        err1(k)=mean(dz.^2);
        err2(k)=mean(abs(dz));
        err3(k)=std(dz);
%         dzerr{k}=ones(size(beads(k).loc.znm))+NaN;
         dzerr{k}=Zint(zglass,zh)-z0f;
   end
end
%
function  f0=getf0Z(loc,p)
frames=loc.frame;
z=loc.znm;
window=round(100/p.dz);
az=abs(z);
nn=loc.phot./(az+p.dz);
[~,maxind]=max(nn);

range=max(1,maxind-window):min(maxind+window,length(z));

zs=double(z(range));
fs=double(frames(range));
fp=fit(zs,fs,'poly1');
f0=fp(0);
end

% function f0=getfsxsy(bead,p)
% roi=15;
% s=size(bead.stack.image);
% dn=(s(1)-roi)/2;
% stackh=bead.stack.image(dn:end-dn,dn:end-dn,:);
% P=callYimingFitter(single(stackh),1,100,4,0,0);
% sx=P(:,5);
% sy=P(:,6);
% phot=P(:,3);
% p.ploton=1;
% figure(77)
% [zas,zn]=stackas2z(sx,sy,bead.stack.framerange',phot,p);
% f0=zas;
% end


function getcoords(a,b,obj)
locsall=obj.locData.getloc({'frame','xnm','ynm','PSFxnm','PSFynm','filenumber','phot'},'position','all','layer',1,'removeFilter','filenumber');
p=obj.getAllParameters;
x=locsall.xnm/p.cam_pixelsize_nm;
y=locsall.ynm/p.cam_pixelsize_nm;
img=myhist2(x,y,1,1,[0 512.5],[0 512.5]);
f1=figure;
imagesc(img');
h=imrect;
pos=wait(h);% [x y wx wy]
delete(f1);
answ=inputdlg({'number of rows','number of columns'},'set tiles', 1,{'2','1'});
if isempty(answ)
    return;
end
nx=str2double(answ{2});
ny=str2double(answ{1});

% pos=pos([2 1 4 3]);
p.Xmin=round(pos(1)); p.Ymin=round(pos(2)); p.Xmax=round(pos(1)+pos(3)); p.Ymax=round(pos(2)+pos(4));
p.Xd=floor(pos(3)/nx);
p.Yd=floor(pos(4)/ny);
obj.setGuiParameters(p)

end


function pard=guidef(obj)
tp=3.6;tmin=4.1;td=4.4;tmax=4.7;
w=0.3;
wp=0.5;
wcb=1.;

pard.beadsource.object=struct('String',{{'RoiManager','Segment'}},'Style','popupmenu','Value',2);
pard.beadsource.position=[1,1];
pard.beadsource.Width=1;

pard.dzt.object=struct('Style','text','String','dz (nm)'); 
pard.dzt.position=[2,1];
pard.dzt.Width=.4;
pard.dz.object=struct('Style','edit','String','10'); 
pard.dz.position=[2,1.4];
pard.dz.Width=.35;


% pard.alignz.object=struct('Style','checkbox','String','Align in z with f0','Value',1); 
% pard.alignz.position=[7,tp];
% pard.alignz.Width=1.2;


pard.zrangeuset.object=struct('String','zrange (nm)','Style','text');
pard.zrangeuset.position=[3,1.5];
pard.zrangeuset.Width=.65;
pard.zrangeuse.object=struct('String','-800 800','Style','edit');
pard.zrangeuse.position=[3,2.15];
pard.zrangeuse.Width=.65;

pard.smoothingt.object=struct('String','Smoothing (frame, zfit)','Style','text');
pard.smoothingt.position=[4,1];
pard.smoothingt.Width=1;
pard.smoothing.object=struct('String','0.02 0.0001','Style','edit');
pard.smoothing.position=[4,2];
pard.smoothing.Width=1;

% pard.spatialcalibration.object=struct('Style','checkbox','String','Spatial calibration','Value',0,'Callback',{{@setvisible,obj}}); 
% pard.spatialcalibration.position=[1,tp];
% pard.spatialcalibration.Width=wcb+.2;
% 
% pard.getcoords.object=struct('String','select','Style','pushbutton','Callback',{{@getcoords,obj}});
% pard.getcoords.position=[1,tp+wcb];
% pard.getcoords.Width=3*w+wp-wcb;
% pard.getcoords.Height=1;

% pard.tt1.object=struct('String','grid','Style','text');
% pard.tt1.position=[2,tp];
% pard.tt1.Width=wp;
% pard.mint.object=struct('String','min','Style','text');
% pard.mint.position=[2,tmin];
% pard.mint.Width=w;
% pard.dxt.object=struct('String','delta','Style','text');
% pard.dxt.position=[2,td];
% pard.dxt.Width=w;
% pard.maxt.object=struct('String','max','Style','text');
% pard.maxt.position=[2,tmax];
% pard.maxt.Width=w;
% 
% pard.Xt.object=struct('String','X (pix)','Style','text');
% pard.Xt.position=[3,tp];
% pard.Xt.Width=wp;
% 
% pard.Xmin.object=struct('String','0','Style','edit');
% pard.Xmin.position=[3,tmin];
% pard.Xmin.Width=w;
% pard.Xd.object=struct('String','512','Style','edit');
% pard.Xd.position=[3,td];
% pard.Xd.Width=w;
% pard.Xmax.object=struct('String','512','Style','edit');
% pard.Xmax.position=[3,tmax];
% pard.Xmax.Width=w;
% 
% pard.Yt.object=struct('String','Y (pix)','Style','text');
% pard.Yt.position=[4,tp];
% pard.Yt.Width=wp;
% 
% pard.Ymin.object=struct('String','0','Style','edit');
% pard.Ymin.position=[4,tmin];
% pard.Ymin.Width=w;
% pard.Yd.object=struct('String','256','Style','edit');
% pard.Yd.position=[4,td];
% pard.Yd.Width=w;
% pard.Ymax.object=struct('String','512','Style','edit');
% pard.Ymax.position=[4,tmax];
% pard.Ymax.Width=w;
% 
% pard.xyoverlapt.object=struct('String','x,y overlap (pix)','Style','text');
% pard.xyoverlapt.position=[5,tp];
% pard.xyoverlapt.Width=1;
% pard.xyoverlap.object=struct('String','10','Style','edit');
% pard.xyoverlap.position=[5,tmax];
% pard.xyoverlap.Width=w;
% 
% pard.zcalc.object=struct('String','z-dependent calibration','Style','checkbox','Value',1,'Callback',{{@setvisible,obj}});
% pard.zcalc.position=[6,tp];
% pard.zcalc.Width=1.5;
% 
% pard.Zt.object=struct('String','Z vals (nm)','Style','text');
% pard.Zt.position=[7,tp];
% pard.Zt.Width=wp;
% pard.Zval.object=struct('String',' 0:1000:3000','Style','edit');
% pard.Zval.position=[7,tmin];
% pard.Zval.Width=5-tp-wp;
% 
% pard.framerangecombinet.object=struct('String','z beyond interval (nm)','Style','text');
% pard.framerangecombinet.position=[8,tp];
% pard.framerangecombinet.Width=1;
% pard.framerangecombine.object=struct('String','100','Style','edit');
% pard.framerangecombine.position=[8,tmax];
% pard.framerangecombine.Width=w;







pard.inputParameters={'cam_pixelsize_nm'};
pard.plugininfo.type='ProcessorPlugin';


end
