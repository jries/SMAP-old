classdef calibrater3DAcomplete<interfaces.DialogProcessor
    properties
        SXY
    end
    methods
        function obj=calibrater3DAcomplete(varargin)   
            obj@interfaces.DialogProcessor(varargin{:}) ;
             obj.showresults=true;
        end
        function initGui(obj)
            beaddistribution_callback(0,0,obj)           
        end
        function out=run(obj,p)
            out=[];
            
            
            locsall=obj.locData.getloc({'frame','xnm','ynm','PSFxnm','PSFynm','filenumber','phot'},'position','all','layer',1,'removeFilter','filenumber');
%             locsall=obj.locData.getloc('frame','xnm','ynm','PSFxnm','PSFynm','filenumber','phot');
            if isempty(locsall.PSFynm)
                error('no PSFy found')
            end
            locsall.PSFxpix=locsall.PSFxnm/p.cam_pixelsize_nm;
            locsall.PSFypix=locsall.PSFynm/p.cam_pixelsize_nm;
            locsall.zfnm=locsall.frame*p.dz;
            
            ax=initaxis(p.resultstabgroup,'found beads');
            htg=uitabgroup(ax.Parent);
                         axbeadss=maketgax(htg,'all');
                         axbeadss.NextPlot='add';
            numfiles=max(locsall.filenumber);
            maxd=p.cam_pixelsize_nm*2;
            locsall.beadnum=zeros(size(locsall.xnm));
            for k=1:numfiles
                ax=maketgax(htg,num2str(k));
                indf=locsall.filenumber==k;
                beadlocs(k)=getBeadLocs(locsall.xnm(indf),locsall.ynm(indf),p);
                [beadnum,numlocs]=associatelocs(beadlocs(k).x,beadlocs(k).y,locsall.xnm(indf),locsall.ynm(indf),maxd);
                beadn=beadnum>0;
                indff=find(indf);
                locsall.beadnum(indff(beadn))=beadnum(beadn)+max(locsall.beadnum);
                plot(axbeadss,beadlocs(k).x,beadlocs(k).y,'o')
            end
            
            %make bead structure
            
%             [ztrue,I0]=getTrueZ(locsall,p);
%              initaxis(p.resultstabgroup,'get true z')
            for k=max(locsall.beadnum):-1:1
                thisbead=(locsall.beadnum==k);
                bead(k).PSFxpix=double(locsall.PSFxpix(thisbead));
                bead(k).PSFypix=double(locsall.PSFypix(thisbead));
                bead(k).phot=double(locsall.phot(thisbead));
                bead(k).xpos=median(double(locsall.xnm(thisbead)));
                bead(k).ypos=median(double(locsall.ynm(thisbead)));
                bead(k).filenumber=double(locsall.filenumber(find(thisbead,1)));
                bead(k).numlocs=sum(thisbead);
                bead(k).frame=double(locsall.frame(thisbead));
                bead(k).zfnm=double(locsall.zfnm(thisbead));
                %derived properties
                bead(k).ztrue=stackas2z(bead(k).PSFxpix,bead(k).PSFypix,bead(k).zfnm,bead(k).phot,0);
%                 bead(k).ztrue
                indz0=find(bead(k).zfnm>bead(k).ztrue,1);
                rangez0=indz0-15:indz0+15;rangez0(rangez0<1)=1;rangez0(rangez0>length(bead(k).zfnm))=length(bead(k).zfnm);
                if isempty(rangez0)
                bead(k).minSx=NaN;
                bead(k).minSy=NaN;
                else
                bead(k).minSx=min(bead(k).PSFxpix(rangez0));
                bead(k).minSy=min(bead(k).PSFypix(rangez0));                    
                end
%                 ztrue(k);
                bead(k).I0=max(bead(k).phot);
            end
            
            % if on glass: correct position based on average, not taking
            % into account position. This might create errors. change
            % later?
            zrangeall=[min(locsall.frame) max(locsall.frame)].*p.dz;
            if p.beaddistribution.Value==1 %glass
                f1ind=find([bead(:).filenumber]==1);
                zf1=robustMean([bead(f1ind).ztrue]);
                if p.ztoframet
                    zpos=p.ztoframe.*p.dz;
                else
                    zpos=zf1;
                end
                zrangeall=zrangeall-zpos;
                for k=1:max([bead(:).filenumber])
                        bi=find([bead(:).filenumber]==k);
                        dzh=robustMean([bead(bi).ztrue])-zf1;
                        for b=1:length(bi)
%                             dzh+zpos
                            bead(bi(b)).ztrue=bead(bi(b)).ztrue-dzh-zpos;
                            bead(bi(b)).zfnm=bead(bi(b)).zfnm-dzh-zpos;
                            bead(bi(b)).zrangeall=zrangeall;
                        end
                end
            end
            
            %sort beads according to X,Y
            
            Xrange=[ p.Xmin:p.Xd:p.Xmax].*p.cam_pixelsize_nm;
            Yrange=[p.Ymin:p.Yd:p.Ymax].*p.cam_pixelsize_nm;
%             Zrange=p.Zmin:p.Zd:p.Zmax;
            Zrange=p.Zval;
            if length(Xrange)==1
                Xrange(2)=p.Xmax*p.cam_pixelsize_nm;
            end
            if length(Yrange)==1
                Yrange(2)=p.Ymax*p.cam_pixelsize_nm;
            end
            if length(Zrange)==1
                Zrange(2)=p.Zmax;
            end
            
            ax=initaxis(p.resultstabgroup,'splines');
            htsplines=uitabgroup(ax.Parent);
            ax=initaxis(p.resultstabgroup,'sx^2-sy^2');
            htsx2sy2=uitabgroup(ax.Parent);
            ax=initaxis(p.resultstabgroup,'z fit');
            htzfit=uitabgroup(ax.Parent);
            ax=initaxis(p.resultstabgroup,'Z(Sx,Sy)');
            hzsx=uitabgroup(ax.Parent);
            axzfits=maketgax(htzfit,'summary');  
            axsxsys=maketgax(htsx2sy2,'summary');  
            axsxs22=maketgax(htsx2sy2,'validation');  
            axsplines=maketgax(htsplines,'summary'); 
            axzlut=maketgax(hzsx,'validation'); 
            
            for X=1:length(Xrange)-1
                for Y=1:length(Yrange)-1
                    indh=[bead(:).xpos]>Xrange(X)&[bead(:).xpos]<Xrange(X+1)&[bead(:).ypos]>Yrange(Y)&[bead(:).ypos]<Yrange(Y+1);
                    if sum(indh)==0
                        continue
                    end
                    bh=cleanupbeads(bead(indh),p);
                    if p.beaddistribution.Value==1
                        Z=1;
                        
                        tgt=[num2str(X) num2str(Y) num2str(Z)];
                        ax=maketgax(htsplines,tgt);                     
                        [SXY(X,Y,Z).spline,indg]=getcleanspline(bh,p);
                        bh=bh(indg);
                        ax=maketgax(htsx2sy2,tgt);
                        SXY(X,Y,Z).Sx2_Sy2=cal_Sx2_Sy2(bh,SXY(X,Y,Z).spline,p);
                        ax=maketgax(htzfit,tgt);
                        SXY(X,Y,Z).fitzpar=cal_fitzpar(bh,p);
                        if p.calculateZSxSy
                            ax=maketgax(hzsx,tgt);
                            SXY(X,Y,Z).splineLUT=cal_splineLUT(SXY(X,Y,Z).spline,p);
                        end
                    else %get gel curves
                        splinez=getcleangel(bead(indh),p);
                        for Z=1:length(splinez)
                            SXY(X,Y,Z).spline=splinez(Z);
                        end
                        
                    end 
                    SXY(X,Y,Z).bead=bh;
                    SXY(X,Y,Z).Xrangeall=Xrange;
                    SXY(X,Y,Z).Yrangeall=Yrange;
                    SXY(X,Y,Z).Zrangeall=Zrange;
                    SXY(X,Y,Z).posind=[X,Y,Z];
                    SXY(X,Y,Z).Xrange=[Xrange(X), Xrange(X+1)];
                    SXY(X,Y,Z).Yrange=[Yrange(Y) ,Yrange(Y+1)];
                    SXY(X,Y,Z).Zrange=[Zrange(Z), Zrange(Z+1)];

                end
            end
            
            
           
            axes(axbeadss)
            sp=SXY(:);
            for k=1:numel(sp)
                if ~isempty(sp(k).bead)
                    xpos=vertcat(sp(k).bead(:).xpos);
                    ypos=vertcat(sp(k).bead(:).ypos);
                plot(xpos,ypos,'k.');hold on;
                end
            end
            for k=1:length(Xrange)
                line([Xrange(k),Xrange(k)],[Yrange(1) Yrange(end)])
            end
            for k=1:length(Yrange)
                line([Xrange(1),Xrange(end)],[Yrange(k) Yrange(k)])
            end
            
            sp=[SXY(:).spline];
%             zt=zrangeall(1):0.01:zrangeall(2);
            z0=zf1-zpos;
            zt=z0+p.zrangeuse(1):0.01:z0+p.zrangeuse(2);

            
            axes(axsplines)
            hold off;for k=1:numel(sp), plot(axsplines,zt,sp(k).x(zt),zt,sp(k).y(zt));hold on;end
            ylim([1 5])
            
            
            sp={SXY(:).Sx2_Sy2};
            s=-30:0.5:30;
            axes(axsxsys)
            hold off;
            for k=1:numel(sp)
                if ~isempty(sp{k})
                plot(axsxsys,sp{k}(s),s);hold on;
                end
            end
            xlabel(axsxsys,'z')
            ylabel(axsxsys,'Sx^2-Sy^2')
            xlim(axsxsys,p.zrangeuse)
            
            sp={SXY(:).fitzpar};
            
            axes(axzfits)
            hold off
            xlabel(axzfits,'z')
            ylabel(axzfits,'Sx, Sy')
            ylim(axzfits,[0 5])
            
            hold off;
            for k=1:numel(sp)
                if ~isempty(sp{k})
                    [sx,sy]=getsxfromzfitpar(zt,sp{k},z0); 
                    plot(axzfits,zt,sx,zt,sy)
                    hold on
                end
            end
            
            
            %cross check and validate with bead positions
%             axes(axzlut);
            
%             hold off
            zt=p.zrangeuse(1):p.dz:p.zrangeuse(2);
            
            for k=1:numel(SXY)
                sxa=vertcat(SXY(k).bead(:).PSFxpix);
                sya=vertcat(SXY(k).bead(:).PSFypix);
                za=vertcat(SXY(k).bead(:).zfnm);
                
                if ~isempty(SXY(k).splineLUT)
                    zha=zfromSXSYLut(SXY(k).splineLUT,sxa,sya);
                    dzba=bindata(za,za-zha,zt,'mean');
                    dzs=bindata(za,za-zha,zt,'std');
                    plot(axzlut,za,za-zha,'.','MarkerSize',2)
                    axzlut.NextPlot='add';
                    plot(axzlut,zt,dzba,'k',zt,dzba+dzs,'k',zt,dzba-dzs,'k')
                end
                    zha2=zfromSx2_Ss2(SXY(k).Sx2_Sy2,sxa,sya);
                     plot(axsxs22,za,za-zha2,'.','MarkerSize',2)
                     axsxs22.NextPlot='add';
                    dzba2=bindata(za,za-zha2,zt,'mean');
                    dzs2=bindata(za,za-zha2,zt,'std');
                    plot(axsxs22,zt,dzba2,'k',zt,dzba2+dzs2,'k',zt,dzba2-dzs2,'k')
            end
            plot(axzlut,p.zrangeuse,[0 0],'k');
            plot(axsxs22,p.zrangeuse,[0 0],'k');
            yrange=[-100 100];
            ylim(axzlut,yrange)
            xlim(axzlut,p.zrangeuse)
            ylim(axsxs22,yrange)
            xlim(axsxs22,p.zrangeuse)
            
%             axes(axsxs22);
%             hold off
%             for k=1:numel(SXY)
%                 for be=1:length(SXY(k).bead)
%                     zh=zfromSx2_Ss2(SXY(k).Sx2_Sy2,SXY(k).bead(be).PSFxpix,SXY(k).bead(be).PSFypix);
%                     
% %                     zh=zfromSXSYLut(SXY(k).splineLUT,SXY(k).bead(be).PSFxpix,SXY(k).bead(be).PSFypix);
%                     zb=SXY(k).bead(be).zfnm;
%                     plot(zb,zh-zb,'.');
%                     hold on
%                 end
%             end
%             ylim([-150 150])
%             xlim(p.zrangeuse)
            
%             za=zfromSXSYLut(SXY(1).splineLUT,locsall.PSFxpix,locsall.PSFypix);
%             zaf=locsall.zfnm;
%             dz=zaf-za;
%             
            
            
            obj.SXY=SXY;
            
            
            %if in gel: make Z.dependent curves
            
            % calculate average curves 
            
            %remove some curves based on intensity, deviation from average
            %Sxy
            
            %recalculate average curves
            
            %calculate LUT
            
%             calculate sx.^2-sy.^2 
            
%             zglass=min(ztrue)
%             if p.check_glassposb
%                 zglass=p.cal_glasspos*p.dz/1000;
%             else
%                 zs=sort(ztrue(~isnan(ztrue)));
%                 zglass=mean(zs(1:3));
%             end
%             ztrue=ztrue-zglass;
%             title(['zglass = ' num2str(zglass)])
%             
%             B0=double(~p.B0);
%             
%             frame=(p.cal_fstart+p.cal_fstop)/2;
%             [sx,sy,zcorr]=getCalibrationCurve(locsall,ztrue,frame,p);
%             initaxis(p.resultstabgroup,'calibration curve fit');
% 
%             fitp=fitCalibrationCurve(sx./p.cam_pixelsize_nm,sy/p.cam_pixelsize_nm,zcorr,frame*p.dz/1000-zglass,B0,p.cal_zrange/2);          
%             obj.outforfit=real(fitp([2 4 5 6 7 8 1 3]));
%             drawnow
% %             obj.cal3D.zglass=zglass;
%             allframes=p.cal_fstart:p.cal_deltaf:p.cal_fstop;
%             
% 
%             initaxis(p.resultstabgroup,'sx2-sy2');
%             obj.outsx2sy2=fitsx2sy2(sx/p.cam_pixelsize_nm,sy/p.cam_pixelsize_nm,zcorr,frame*p.dz/1000-zglass,p.cal_zrange/2);
%             
% %             recgui.initaxis(p.resultstabgroup,'sx-sy');
% %             plot(zcorr,sx-sy,'.')
%         
%             for k=1:length(allframes)
%                     initaxis(p.resultstabgroup,['c' num2str(k)]);
%                 obj.cal3D(k).zglass=zglass;
%                 frame=allframes(k);
%                 [sx,sy,zcorr]=getCalibrationCurve(locsall,ztrue,frame,p);
% %                 recgui.initaxis(p.resultstabgroup,'calibration curve fit');
%                 fitp=fitCalibrationCurve(sx./p.cam_pixelsize_nm,sy/p.cam_pixelsize_nm,zcorr,frame*p.dz/1000-zglass,B0,p.cal_zrange/2);          
%                 obj.cal3D(k).fit=real(fitp([2 4 5 6 7 8 1 3]));
%                 obj.cal3D(k).midpoint=-fitp(9);
%                 obj.cal3D(k).sx2sy2=fitsx2sy2(sx./p.cam_pixelsize_nm,sy./p.cam_pixelsize_nm,zcorr,frame*p.dz/1000-zglass,p.cal_zrange/2);
%                 
%                 drawnow
%                 
%             end
            
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function save_SXY(a,b,obj)
fileout=obj.locData.files.file(1).name;
fn=[fileout(1:end-4) '_3DAcal.mat'];
[f,p]=uiputfile(fn);
if f
    SXY=obj.SXY;
    save([p f],'SXY')      
    obj.setPar('cal3D_file',[p f]);
end


end

function fitzpar=cal_fitzpar(b,p)
 z=vertcat(b(:).zfnm);
Sx=vertcat(b(:).PSFxpix);
Sy=vertcat(b(:).PSFypix);
% ztrue=robustMean([b(:).ztrue]);
% zrange=spline.maxmaxrange;
% zrange=zrange+[300 -300];

midpoint=robustMean([b(:).ztrue]);
zrange=[midpoint+p.fitzrange(1) midpoint+p.fitzrange(2)];
ax=gca;hold off;
fitzpar=getzfitpar(Sx,Sy,z,zrange,midpoint,~p.fitzB0,ax);

end

function splineLUT=cal_splineLUT(spline,p)
srange=p.Smin:p.Sd:p.Smax;
splineLUT=gets2z(spline,srange);
% splineLUT=[];
end

function ax=maketgax(htg,txt)
ht=uitab(htg,'Title',txt);
ax=axes('Parent',ht);
end

function sxp=cal_Sx2_Sy2(b,spline,p)


 z=vertcat(b(:).zfnm);
Sx=vertcat(b(:).PSFxpix);
Sy=vertcat(b(:).PSFypix);
% figure(88)
hold off
sxp=fitsx2sy2(Sx,Sy,z,spline.maxmaxrange,gca);
end

function bo=cleanupbeads(b,p)
%filter based on Sx, Sy and phot
minSx=[b(:).minSx];
minSy=[b(:).minSy];
% zfnm=[b(:).zfnm];
[mSx,smSx]=robustMean(minSx);
[mSy,smSy]=robustMean(minSy);
phot=[b(:).I0];
[mp,sp]=robustMean(phot);
ztrue=[b(:).ztrue];
indgood=minSx<mSx+2*smSx & minSy<mSy+2*smSy & phot<mp+2*sp&~isnan(minSx)&~isnan(minSy)&imag(ztrue)==0;
indgood=indgood&ztrue>-10000 &ztrue<30000;

if sum(indgood)==0
    indgood=true(size(minSx));
end
% correctindividual=0; %there seems to be no tilt or similar
% 
% if correctindividual
%     
%     for k=1:length(b)
%         b(k).zfnm=b(k).zfnm-b(k).ztrue+zm;
%     end
% end

bo=b(indgood);
end
function [s,indgood2]=getcleanspline(b,p)

zm=robustMean([b(:).ztrue]);


 z=vertcat(b(:).zfnm);
Sx=vertcat(b(:).PSFxpix);
Sy=vertcat(b(:).PSFypix);
indz=z>zm+p.zrangeuse(1)&z<zm+p.zrangeuse(2);
z=z(indz);
Sx=Sx(indz);
Sy=Sy(indz);

splinex=getspline(Sx,z,1./Sx,p.splinepar);
spliney=getspline(Sy,z,1./Sy,p.splinepar);

% figure(88);plot(z,Sx,'.')
%how good does it fit with spline?
% figure(99)
hold off
for k=length(b):-1:1
    w=(b(k).phot).^2;
%     w=1;
    errh=(b(k).PSFxpix-splinex(b(k).zfnm)).^2.*w+(b(k).PSFypix-spliney(b(k).zfnm)).^2.*w;
    err(k)=sqrt(sum(errh)/sum(w));
%     plot(b(k).zfnm,b(k).PSFxpix,'y.')
%     hold on
%     plot(b(k).zfnm,b(k).PSFypix,'y.')
%     plot(b(k).zfnm,w/max(w))
    
end

%  
[em,es]=robustMean(err);
if isnan(es)
    es=0;
end
indgood2=err<em+2*es;
if sum(indgood2)==0
    indgood2=true(length(b),1);
end
zg=vertcat(b(indgood2).zfnm);
indz=zg>zm+p.zrangeuse(1)&zg<zm+p.zrangeuse(2);
zg=zg(indz);
sxg=vertcat(b(indgood2).PSFxpix);
syg=vertcat(b(indgood2).PSFypix);
sxg=sxg(indz);
syg=syg(indz);

splinex2=getspline(sxg,zg,1./(sxg-splinex(zg)).^2,p.splinepar);
spliney2=getspline(syg,zg,1./(syg-spliney(zg)).^2,p.splinepar);

% zt=b(1).zrangeall(1):0.01:b(1).zrangeall(2);
zt=min(zg):0.01:max(zg);

for k=1:length(b)
    if indgood2(k)
        t='r.';
    else
        t='g.';
    end
    plot(b(k).zfnm,b(k).PSFxpix,t)
    hold on
    plot(b(k).zfnm,b(k).PSFypix,t)
    
end

plot(zt,splinex(zt),'b')
plot(zt,spliney(zt),'b')
plot(zt,splinex2(zt),'k')
plot(zt,spliney2(zt),'k')

xlim([zt(1) zt(end)])
drawnow
s.x=splinex2;s.y=spliney2;
s.zrange=[zt(1) zt(end)];


zr=s.zrange(1):0.01:s.zrange(2);
midp=round(length(zr)/2);

[~,ind1x]=max(s.x(zr(1:midp)));
[~,ind1y]=max(s.y(zr(1:midp)));
[~,ind2x]=max(s.x(zr(midp:end)));
[~,ind2y]=max(s.y(zr(midp:end)));

z1=max(zr(ind1x),zr(ind1y));
z2=min(zr(ind2x+midp-1),zr(ind2y+midp-1));

s.maxmaxrange=[z1 z2];
end
function spline=getspline(S,z,w,p)
if nargin<4
    p=0.96;
end
% 
[zs,zind]=sort(z);Ss=S(zind);ws=w(zind);
spline=fit(zs,Ss,'smoothingspline','Weights',ws,'Normalize','on','SmoothingParam',p); 
% spline=fit(zs,Ss,'smoothingspline','Weights',ws,'Normalize','on'); 
% spline=fit(zs,Ss,'smoothingspline','Weights',ws,'SmoothingParam',0.99); %weigh smaller more
% [spline,a,b]=fit(zs,Ss,'smoothingspline','Weights',ws); %weigh smaller more
% ds=Ss-spline(zs);
% spline2=fit(zs,Ss,'smoothingspline','Weights',1./ds.^2); %weigh smaller more
end

% 
% function savecalfile(a,b,obj)
% fn=[obj.locData.files.file(1).name(1:end-4) '_3DAcal.mat'];
% [f,p]=uiputfile(fn);
% if f
%     outforfit=obj.outforfit;
%     cal3D=obj.cal3D;
%     outsx2sy2=obj.outsx2sy2;
%     save([p f],'outforfit','cal3D','outsx2sy2')
% end
% end
% 
% function fitpsx=fitsx2sy2(sx,sy,zcorr,framez,zrange)
% 
% indf=abs(zcorr-framez)<zrange;
% if sum(indf)>5
% fitpsx=fit(sx(indf).^2-sy(indf).^2,zcorr(indf),'poly3','Robust','LAR');
% 
% % fitpsx=polyfit(zcorr,sx.^2-sy.^2,4);
% plot(sx.^2-sy.^2,zcorr,'r.')
% hold on
% plot(sx(indf).^2-sy(indf).^2,zcorr(indf),'b.')
% 
% sxsort=sort(sx.^2-sy.^2);
% zsort=feval(fitpsx,sxsort);
% 
% plot(sxsort,zsort,'k')
% % plot(zcorr,polyval(fitpsx,zcorr),'.')
% hold off
% xlim([-6 6])
% else
%     fitpsx=zeros(2,1);
% end
% end

% function fitp=fitCalibrationCurve(sx,sy,z,framez,B0,zrange,startp)
% if isempty(sx)
%     fitp=zeros(9,1);
% else
% if nargin<7
%     startp=[    0.3    1.0    1.0000  0   0        0         0  0.307   -framez];
% end
% indf=abs(z-framez)<zrange;
% fitp=lsqnonlin(@sbothfromsigmaerr,startp,[],[],[],[z(indf) z(indf)],[sx(indf) sy(indf)],B0);
% fitp=real(fitp);
% 
% subplot('Position',[0.05,0.65,.9,.3])
% plot(z,sx,'c.',z,sy,'m.')
% hold on;
% plot(z(indf),sx(indf),'b.',z(indf),sy(indf),'r.')
% sxf=sigmafromz(fitp([1 2 4 6 8 9]),z,B0);
% plot(z,sxf,'k.')
% fpy=fitp([1 3 5 7 8 9]);
% fpy(5)=-fpy(5);
% syf=sigmafromz(fpy,z,B0);
% plot(z,syf,'k.')
% hold off
% title(num2str(fitp(:)',2))
% ylim([0 max(max(sxf),max(syf))])
% 
% subplot('Position',[0.05,0.45,.9,.15])
% plot(z(indf),sx(indf)-sxf(indf),'b.')
% hold on
% plot(z(indf),0*sy(indf),'k.')
% hold off
% ylim([-1 1]*.7)
% 
% subplot('Position',[0.05,0.25,.9,.15])
% plot(z(indf),sy(indf)-syf(indf),'r.')
% hold on
% plot(z(indf),0*sy(indf),'k.')
% hold off
% ylim([-1 1]*.7)
% 
% subplot('Position',[0.05,0.05,.9,.15])
% end
% end
function [sx,sy,zcorr]=getCalibrationCurve(locs,ztrue,frame,p)
frameregion=p.cal_framerange/2;        

df=locs.frame-frame;

indf=abs(df)<=frameregion;
indbead=locs.beadnum>0;

sx=double(locs.PSFxnm(indf&indbead));
sy=double(locs.PSFynm(indf&indbead));
fn=locs.frame(indf&indbead);
beadnum=locs.beadnum(indf&indbead);
dff=df(indf&indbead);

zbead=ztrue(beadnum);
zcorr=zbead-dff*p.dz;
% 
%  recgui.initaxis(p.resultstabgroup,'calibration curve')
% plot(zcorr,sx,'.',zcorr,sy,'.')

indg=~isnan(zcorr);
zcorr=zcorr(indg);sx=sx(indg);sy=sy(indg);
end

% function [z,I0]=getTrueZ(locs,p)
% mmax=max(locs.beadnum);
% z=zeros(mmax,1);
% I0=zeros(mmax,1);
%  initaxis(p.resultstabgroup,'get true z')
% 
% for k=1:mmax
%     ind=locs.beadnum==k;
% %     zf=locs.frame(ind)*p.dz/1000;
%     [zas,zn]=stackas2z(locs.PSFxpix(ind),locs.PSFypix(ind),locs.zfnm(ind),locs.phot(ind),p.showresults);
%     z(k)=zas;
%     I0(k)=max(locs.phot(ind));
% %     waitforbuttonpress
% end
%  w=warning('off');
%  warning(w);
% end

function calibrateAstig3D(locs,p)
global zt sxf syf
sxt=double(locs.PSFxnm)/p.cam_pixelsize_nm;syt=double(locs.PSFynm)/p.cam_pixelsize_nm;framet=double(locs.frame);
zt=framet*p.dz;

B0=double(~p.B0);


recgui.initaxis(p.resultstabgroup,'select range');
plot(framet,sxt,'ro')
hold on
plot(framet,syt,'bo')
hold off
title('select range (two points)')
[indi,y]=ginput(2);

% indi=round(indum/p.dz*1000);

%rangef=round(max(1,min(indi))):round(min(max(indi),length(sxt)));
rangef=round(max(min(indi))):round(min(max(indi)));

range=find(framet>=rangef(1),1,'first'):find(framet<=rangef(end),1,'last');

sx=sxt(range);sy=syt(range);frame=framet(range);
z=frame*p.dz;



[~,ix]=min(sx);
[~,iy]=min(sy);
% 
midpreal=(z(ix)+z(iy))/2;


indframez0=find(p.framez0<=framet,1,'first');
midp=zt(indframez0);

z=z-midp;zt=zt-midp;
% recgui.initaxis(p.resultstabgroup,'fitted points');

plot(zt,sxt,'r.')
hold on
plot(zt,syt,'b.')
plot(z,sx,'ro')
plot(z,sy,'bo')

mpf=midpreal-midp;
% parx= [d sx0 sy0 Ax Ay Bx By g mp]
startp=[    0.3    1.0    1.0000  0   0        0         0  0.307   -mpf];

% fitp=lsqcurvefit(@sbothfromsigma,startp,[z z],[sx sy])
fitp=lsqnonlin(@sbothfromsigmaerr,startp,[],[],[],[z z],[sx sy],B0);

% fitfunction=@(par,x)sbothfromsigmafit(par,x);
% fit2=fit([z; z],[sx; sy],fitfunction);

sxf=sigmafromz(fitp([1 2 4 6 8 9]),zt,B0);
plot(zt,sxf,'k')
fpy=fitp([1 3 5 7 8 9]);
fpy(5)=-fpy(5);
syf=sigmafromz(fpy,zt,B0);
plot(zt,syf,'k')
hold off

axis tight

% 
% figure(2)
% plot(sx,sy,'.')
% hold on
% plot(sxf,syf,'r')
% plot(syf,sxf,'c')
% hold off


%zpar=[sigma0x,Ax,Ay,Bx,By,gamma,d,sigma0y)
outforfit=real(fitp([2 4 5 6 7 8 1 3]));
button = questdlg('Is the fit good?','3D astigmatism bead calibration','Refit','Save','Cancel','Refit') ;
switch button
    case 'Refit'
        calibrateAstig3D(locs,p)
    case 'Save'
        fn=[p.fileout(1:end-4) '_3DAcal.mat']
        [f,p]=uiputfile(fn);
        if f
            save([p f],'outforfit')
        end
end
end
function s=sigmafromz(par,z,B0)
global gamma
par=real(par);
% parx= [d sx0 Ax Bx g mp]
s0=par(2);d=par(1);A=par(3);B=par(4)*B0;g=par(5);mp=par(6);

% s=s0*sqrt(1+(z-g+mp).^2/d^2);
s=s0*sqrt(1+(z-g+mp).^2/d^2+A*(z-g+mp).^3/d^3+B*(z-g+mp).^4/d^4);
s=real(s);
end


function s=sbothfromsigma(par,z,B0)
% parx= [d sx0 sy0 Ax Ay Bx By g mp]
par=real(par);
px=par([1 2 4 6 8 9]);
py=par([1 3 5 7 8 9]);
 py(5)=-py(5);
zh=z(:,1);
s=real([sigmafromz(px,zh,B0) sigmafromz(py,zh,B0)]);
end

function err=sbothfromsigmaerr(par,z,sx,B0)
sf=sbothfromsigma(par,z,B0);
err=sf-sx;

% stderr=std(err);
err=err./sqrt(abs(err));
end

function beaddistribution_callback(a,b,obj)
gel={'Zt','Zval','framerangecombinet','framerangecombine'};
glass={'ztoframe','ztoframet'};
p=obj.getSingleGuiParameter('beaddistribution');
switch p.Value 
    case 1%glass
        tgel='off';tglass='on';
    case 2 %gel
        tgel='on';tglass='off';
end
for k=1:length(gel)
    obj.guihandles.(gel{k}).Visible=tgel;
end
for k=1:length(glass)
    obj.guihandles.(glass{k}).Visible=tglass;
end
end


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
answ=inputdlg({'number of rows','number of columns'},'set tiles', 1,{'3','3'});
if isempty(answ)
    return;
end
nx=str2double(answ{1});
ny=str2double(answ{2});

% pos=pos([2 1 4 3]);
p.Xmin=round(pos(1)); p.Ymin=round(pos(2)); p.Xmax=round(pos(1)+pos(3)); p.Ymax=round(pos(2)+pos(4));
p.Xd=floor(pos(3)/nx);
p.Yd=floor(pos(4)/ny);
obj.setGuiParameters(p)

end

function pard=guidef(obj)
pard.dzt.object=struct('Style','text','String','dz (nm)'); 
pard.dzt.position=[1,1];
pard.dzt.Width=.5;
pard.dz.object=struct('Style','edit','String','50'); 
pard.dz.position=[1,1.5];
pard.dz.Width=.5;

pard.beaddistribution.object=struct('String',{{'Beads on Glass','Beads in Gel'}},'Style','popupmenu','Callback',{{@beaddistribution_callback,obj}});
pard.beaddistribution.position=[2,1];
pard.beaddistribution.Width=2;

tp=3.1;tmin=3.6;td=3.9;tmax=4.2;
w=0.3;
wp=0.5;
pard.tt1.object=struct('String','grid','Style','text');
pard.tt1.position=[2,tp];
pard.tt1.Width=wp;
pard.mint.object=struct('String','min','Style','text');
pard.mint.position=[2,tmin];
pard.mint.Width=w;
pard.dxt.object=struct('String','delta','Style','text');
pard.dxt.position=[2,td];
pard.dxt.Width=w;
pard.maxt.object=struct('String','max','Style','text');
pard.maxt.position=[2,tmax];
pard.maxt.Width=w;

pard.Xt.object=struct('String','X (pix)','Style','text');
pard.Xt.position=[3,tp];
pard.Xt.Width=wp;

pard.Xmin.object=struct('String','0','Style','edit');
pard.Xmin.position=[3,tmin];
pard.Xmin.Width=w;
pard.Xd.object=struct('String','128','Style','edit');
pard.Xd.position=[3,td];
pard.Xd.Width=w;
pard.Xmax.object=struct('String','512','Style','edit');
pard.Xmax.position=[3,tmax];
pard.Xmax.Width=w;

pard.Yt.object=struct('String','Y (pix)','Style','text');
pard.Yt.position=[4,tp];
pard.Yt.Width=wp;

pard.Ymin.object=struct('String','0','Style','edit');
pard.Ymin.position=[4,tmin];
pard.Ymin.Width=w;
pard.Yd.object=struct('String','128','Style','edit');
pard.Yd.position=[4,td];
pard.Yd.Width=w;
pard.Ymax.object=struct('String','512','Style','edit');
pard.Ymax.position=[4,tmax];
pard.Ymax.Width=w;

pard.Zt.object=struct('String','Z vals (nm)','Style','text');
pard.Zt.position=[5,tp];
pard.Zt.Width=wp;
pard.Zval.object=struct('String',' 0:1000:3000','Style','edit');
pard.Zval.position=[5,tmin];
pard.Zval.Width=5-tp-wp;

% pard.Zmin.object=struct('String','0','Style','edit');
% pard.Zmin.position=[5,tmin];
% pard.Zmin.Width=w;
% pard.Zd.object=struct('String','1000','Style','edit');
% pard.Zd.position=[5,td];
% pard.Zd.Width=w;
% pard.Zmax.object=struct('String','3000','Style','edit');
% pard.Zmax.position=[5,tmax];
% pard.Zmax.Width=w;

pard.splinet.object=struct('String','Spline','Style','text');
pard.splinet.position=[6,1];
pard.splinet.Width=1;
pard.splinepart.object=struct('String','smoothing par','Style','text');
pard.splinepart.position=[6,2];
pard.splinepart.Width=1;
pard.splinepar.object=struct('String','0.95','Style','edit');
pard.splinepar.position=[6,3];
pard.splinepar.Width=1;

pard.fitzt.object=struct('String','fit z MLE:','Style','text');
pard.fitzt.position=[7,1];
pard.fitzt.Width=1;
pard.fitzranget.object=struct('String','z range (nm)','Style','text');
pard.fitzranget.position=[7,2];
pard.fitzranget.Width=1;
pard.fitzrange.object=struct('String','-500 500','Style','edit');
pard.fitzrange.position=[7,3];
pard.fitzrange.Width=1;
pard.fitzB0.object=struct('String','B0=0','Style','checkbox','Value',0);
pard.fitzB0.position=[7,4];
pard.fitzB0.Width=1;

pard.calculateZSxSy.object=struct('String','Z(Sx,Sy,X,Y,Z)','Style','checkbox','Value',0);
pard.calculateZSxSy.position=[8,1];
pard.calculateZSxSy.Width=1.;

pard.St.object=struct('String','S (pix)','Style','text');
pard.St.position=[8,tp];
pard.St.Width=wp;
pard.Smin.object=struct('String','0','Style','edit');
pard.Smin.position=[8,tmin];
pard.Smin.Width=w;
pard.Sd.object=struct('String','.03','Style','edit');
pard.Sd.position=[8,td];
pard.Sd.Width=w;
pard.Smax.object=struct('String','4','Style','edit');
pard.Smax.position=[8,tmax];
pard.Smax.Width=w;

pard.framerangecombinet.object=struct('String','frames to combine','Style','text');
pard.framerangecombinet.position=[3,1];
pard.framerangecombinet.Width=1;
pard.framerangecombine.object=struct('String','10','Style','edit');
pard.framerangecombine.position=[3,2];
pard.framerangecombine.Width=.5;

pard.ztoframet.object=struct('String','set z to frame','Style','checkbox');
pard.ztoframet.position=[3,1];
pard.ztoframet.Width=1;
pard.ztoframe.object=struct('String','21','Style','edit');
pard.ztoframe.position=[3,2];
pard.ztoframe.Width=.5;


pard.zrangeuset.object=struct('String','zrange to use (nm)','Style','text');
pard.zrangeuset.position=[4,1];
pard.zrangeuset.Width=1;
pard.zrangeuse.object=struct('String','-800 800','Style','edit');
pard.zrangeuse.position=[4,2];
pard.zrangeuse.Width=1;

pard.getcoords.object=struct('String','select','Style','pushbutton','Callback',{{@getcoords,obj}});
pard.getcoords.position=[4,4.5];
pard.getcoords.Width=.5;
pard.getcoords.Height=2;

pard.save.object=struct('String','save','Style','pushbutton','Callback',{{@save_SXY,obj}});
pard.save.position=[9,1];
pard.save.Width=1;
pard.save.Height=1;
% pard.savebutton.object=struct('String','save','Style','pushbutton');
% pard.savebutton.position=[7,3];
pard.inputParameters={'cam_pixelsize_nm'};
pard.plugininfo.type='ProcessorPlugin';


end
