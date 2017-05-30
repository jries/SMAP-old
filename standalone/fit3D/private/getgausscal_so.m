function [gauss,indgood]=getgausscal_so(beads,p)
    curves=[];
    for B=length(beads):-1:1
        beadz0=(beads(B).f0)*p.dz;
        
        if contains(p.zcorr,'astig')||contains(p.zcorr,'corr')
             beadz=(beads(B).loc.frames*p.dz)-beadz0;
        else 
            beadz=(beads(B).loc.frames-p.midpoint)*p.dz;
        end
           sx=beads(B).loc.PSFxpix; 
           sy=beads(B).loc.PSFypix; 
           z=beadz; phot=beads(B).loc.phot;  
        inzr=z>=p.gaussrange(1)&z<=p.gaussrange(2);
        curves(B).sx=double(sx(inzr));
        curves(B).sy=double(sy(inzr));
        curves(B).z=double(z(inzr));
        curves(B).phot=double(phot(inzr));
        curves(B).xpos=beads(B).pos(1);
        curves(B).ypos=beads(B).pos(2);

    end
    
    %get calibrations
    bh=curves;
%     tgt=[num2str(X) num2str(Y) num2str(Z) num2str(p.iter)];
    ax=axes(uitab(p.tabgroup,'Title','sx,sy'));
     [spline,indg]=getcleanspline(curves,p);
%      indgoodc(indgoodc)=indgoodc(indgoodc)&indg;
    indgood=indg;
    bh=bh(indg);
    
%     SXY(Z)=getsxyinit(p,X,Y,Z);
     gauss.spline=spline;
    
        ax=axes(uitab(p.tabgroup,'Title','sx^2-sy^2'));
        gauss.Sx2_Sy2=cal_Sx2_Sy2(bh,p);
   drawnow
    
        ax=axes(uitab(p.tabgroup,'Title','Gauss: fit z'));
        gauss.fitzpar=cal_fitzpar(bh,p);
    drawnow
end


function [s,indgood2]=getcleanspline(curves,p)

% zm=robustMean([b(:).ztrue]);


 za=vertcat(curves(:).z);
Sxa=vertcat(curves(:).sx);
Sya=vertcat(curves(:).sy);
% indz=za>zm+p.zrangeuse(1)&za<zm+p.zrangeuse(2);
indz=za>p.gaussrange(1)&za<p.gaussrange(2);
z=za(indz);
Sx=Sxa(indz);
Sy=Sya(indz);

% Sxs=smoothn({z,Sx});
warning('off','curvefit:fit:iterationLimitReached');
splinex=fit(z,Sx,'poly6','Robust','LAR','Normalize','on');
spliney=fit(z,Sy,'poly6','Robust','LAR','Normalize','on');
warning('on','curvefit:fit:iterationLimitReached');
hold off
for k=length(curves):-1:1
    w=(curves(k).phot);
%     w=1;
    zh=curves(k).z;
    indzh=zh>p.gaussrange(1)&zh<p.gaussrange(2);
    w=w.*indzh;
    errh=(curves(k).sx-splinex(zh)).^2.*w+(curves(k).sy-spliney(zh)).^2.*w;
    err(k)=sqrt(sum(errh)/sum(w));
    errh2=(curves(k).sx-splinex(zh)).*w./curves(k).sx;
    errh3=(curves(k).sy-spliney(zh)).*w./curves(k).sy;
    err2(k)=abs(sum(errh2)/sum(w));
    err3(k)=abs(sum(errh3)/sum(w));
%     err(k)=sqrt(robustMean(errh));
%     plot(b(k).zfnm,b(k).PSFxpix,'y.')
%     hold on
%     plot(b(k).zfnm,b(k).PSFypix,'y.')
%     plot(b(k).zfnm,w/max(w))
    
end

%  
[em,es]=robustMean(err);
% [em2,es2]=robustMean(err2);
% [em3,es3]=robustMean(err3);
if isnan(es), es=em; end
% if isnan(es2), es2=em2; end
% if isnan(es3), es3=em3; end
indgood2=err<em+2*es & err2+err3<.15;
if sum(indgood2)==0
    indgood2=true(1,length(curves));
end
zg=vertcat(curves(indgood2).z);
indz=zg>p.gaussrange(1)&zg<p.gaussrange(2);
zg=zg(indz);
sxg=vertcat(curves(indgood2).sx);
syg=vertcat(curves(indgood2).sy);
sxg=sxg(indz);
syg=syg(indz);

splinex2=getspline(sxg,zg,1./(abs(sxg-splinex(zg))+.1));
spliney2=getspline(syg,zg,1./(abs(syg-spliney(zg))+.1));
% splinex2=getspline(sxg,zg,1,p.splinepar);
% spliney2=getspline(syg,zg,1,p.splinepar);
% zt=b(1).zrangeall(1):0.01:b(1).zrangeall(2);
zt=min(zg):0.01:max(zg);

for k=1:length(curves)
    if indgood2(k)
        t='r.';
    else
        t='g.';
    end
    plot(curves(k).z,curves(k).sx,t)
    hold on
    plot(curves(k).z,curves(k).sy,t)
    
end

plot(zt,splinex(zt),'b')
plot(zt,spliney(zt),'b')
plot(zt,splinex2(zt),'k')
plot(zt,spliney2(zt),'k')

xlim([zt(1) zt(end)])
ylim([0 min(5,max(max(splinex2(zt)),max(spliney2(zt))))])
drawnow
s.x=splinex2;s.y=spliney2;
s.zrange=[zt(1) zt(end)];


zr=s.zrange(1):1:s.zrange(2);
midp=round(length(zr)/8);


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
        fn=[p.fileout(1:end-4) '_3DAcal.mat'];
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


function sxp=cal_Sx2_Sy2(b,p)


z=vertcat(b(:).z);
Sx=vertcat(b(:).sx);
Sy=vertcat(b(:).sy);
% figure(88)
hold off
sxp=fitsx2sy2_so(Sx,Sy,z,p.gaussrange,gca);
% sxp.ztruepos=p.ztruepos;
xlim([p.gaussrange])
end

function fitzpar=cal_fitzpar(b,p)
 z=vertcat(b(:).z);
Sx=vertcat(b(:).sx);
Sy=vertcat(b(:).sy);
% ztrue=robustMean([b(:).ztrue]);
% zrange=spline.maxmaxrange;
% zrange=zrange+[300 -300];

% midpoint=robustMean([b(:).ztrue]);
% ztruepos=0;
zrange=[p.gaussrange(1) p.gaussrange(2)];
ax=gca;hold off;
fitzpar=getzfitpar(Sx,Sy,z,zrange,0,true,ax);

end

function fitpsx=fitsx2sy2_so(sx,sy,z,zrange,ax)

% indf=abs(z)<zrange;
indf=z>zrange(1)&z<zrange(2);

ds=sx.^2-sy.^2;
q=quantile(ds,[0.01 0.99]);
inds=ds>q(1)&ds<q(2);
indf=indf&inds;


if sum(indf)>5
% fitpsx=fit(sx(indf).^2-sy(indf).^2,z(indf),'smoothingSpline','Normalize','on','SmoothingParam',0.95);
warning('off','curvefit:fit:iterationLimitReached');
fitpsx=fit(sx(indf).^2-sy(indf).^2,z(indf),'poly6','Robust','LAR','Normalize','on');
% lastwarn
% yyaxis right
% hold on
% % fitpsx=polyfit(zcorr,sx.^2-sy.^2,4);
% plot(z,sx.^2-sy.^2,'r.')
if nargin>4
    
    plot(ax,z(indf),sx(indf).^2-sy(indf).^2,'.')
    hold on
    sxsort=sort(sx.^2-sy.^2);
    zsort=feval(fitpsx,sxsort);

    plot(ax,zsort,sxsort,'k')
    ylabel(ax,'sx^2-sy^2')
    % plot(zcorr,polyval(fitpsx,zcorr),'.')

end
% xlim([-6 6])
else
    fitpsx=zeros(2,1);
end
end
