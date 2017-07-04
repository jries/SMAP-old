function [gauss]=getgausscal_so(curves,p)

    
%         p.ax=axes(uitab(p.tabgroup,'Title','sx^2-sy^2'));
        p.ax=p.ax_sxsy;
        gauss.Sx2_Sy2=cal_Sx2_Sy2(curves,p);
   
        p.ax=p.ax_z;
%         p.ax=axes(uitab(p.tabgroup,'Title','Gauss: fit z'));
        gauss.fitzpar=cal_fitzpar(curves,p);
    drawnow
end



% function calibrateAstig3D(locs,p)
% global zt sxf syf
% sxt=double(locs.PSFxnm)/p.cam_pixelsize_nm;syt=double(locs.PSFynm)/p.cam_pixelsize_nm;framet=double(locs.frame);
% zt=framet*p.dz;
% 
% B0=double(~p.B0);
% 
% 
% recgui.initaxis(p.resultstabgroup,'select range');
% plot(framet,sxt,'ro')
% hold on
% plot(framet,syt,'bo')
% hold off
% title('select range (two points)')
% [indi,y]=ginput(2);
% 
% % indi=round(indum/p.dz*1000);
% 
% %rangef=round(max(1,min(indi))):round(min(max(indi),length(sxt)));
% rangef=round(max(min(indi))):round(min(max(indi)));
% 
% range=find(framet>=rangef(1),1,'first'):find(framet<=rangef(end),1,'last');
% 
% sx=sxt(range);sy=syt(range);frame=framet(range);
% z=frame*p.dz;
% 
% 
% 
% [~,ix]=min(sx);
% [~,iy]=min(sy);
% % 
% midpreal=(z(ix)+z(iy))/2;
% 
% 
% indframez0=find(p.framez0<=framet,1,'first');
% midp=zt(indframez0);
% 
% z=z-midp;zt=zt-midp;
% % recgui.initaxis(p.resultstabgroup,'fitted points');
% 
% plot(zt,sxt,'r.')
% hold on
% plot(zt,syt,'b.')
% plot(z,sx,'ro')
% plot(z,sy,'bo')
% 
% mpf=midpreal-midp;
% % parx= [d sx0 sy0 Ax Ay Bx By g mp]
% startp=[    0.3    1.0    1.0000  0   0        0         0  0.307   -mpf];
% 
% % fitp=lsqcurvefit(@sbothfromsigma,startp,[z z],[sx sy])
% fitp=lsqnonlin(@sbothfromsigmaerr,startp,[],[],[],[z z],[sx sy],B0);
% 
% % fitfunction=@(par,x)sbothfromsigmafit(par,x);
% % fit2=fit([z; z],[sx; sy],fitfunction);
% 
% sxf=sigmafromz(fitp([1 2 4 6 8 9]),zt,B0);
% plot(zt,sxf,'k')
% fpy=fitp([1 3 5 7 8 9]);
% fpy(5)=-fpy(5);
% syf=sigmafromz(fpy,zt,B0);
% plot(zt,syf,'k')
% hold off
% 
% axis tight
% 
% % 
% % figure(2)
% % plot(sx,sy,'.')
% % hold on
% % plot(sxf,syf,'r')
% % plot(syf,sxf,'c')
% % hold off
% 
% 
% %zpar=[sigma0x,Ax,Ay,Bx,By,gamma,d,sigma0y)
% outforfit=real(fitp([2 4 5 6 7 8 1 3]));
% button = questdlg('Is the fit good?','3D astigmatism bead calibration','Refit','Save','Cancel','Refit') ;
% switch button
%     case 'Refit'
%         calibrateAstig3D(locs,p)
%     case 'Save'
%         fn=[p.fileout(1:end-4) '_3DAcal.mat'];
%         [f,p]=uiputfile(fn);
%         if f
%             save([p f],'outforfit')
%         end
% end
% end

% function s=sigmafromz(par,z,B0)
% global gamma
% par=real(par);
% % parx= [d sx0 Ax Bx g mp]
% s0=par(2);d=par(1);A=par(3);B=par(4)*B0;g=par(5);mp=par(6);
% 
% % s=s0*sqrt(1+(z-g+mp).^2/d^2);
% s=s0*sqrt(1+(z-g+mp).^2/d^2+A*(z-g+mp).^3/d^3+B*(z-g+mp).^4/d^4);
% s=real(s);
% end

% function s=sbothfromsigma(par,z,B0)
% % parx= [d sx0 sy0 Ax Ay Bx By g mp]
% par=real(par);
% px=par([1 2 4 6 8 9]);
% py=par([1 3 5 7 8 9]);
%  py(5)=-py(5);
% zh=z(:,1);
% s=real([sigmafromz(px,zh,B0) sigmafromz(py,zh,B0)]);
% end

% function err=sbothfromsigmaerr(par,z,sx,B0)
% sf=sbothfromsigma(par,z,B0);
% err=sf-sx;
% 
% % stderr=std(err);
% err=err./sqrt(abs(err));
% end


function sxp=cal_Sx2_Sy2(b,p)


z=vertcat(b(:).z);
Sx=vertcat(b(:).sx);
Sy=vertcat(b(:).sy);
% figure(88)
hold off
sxp=fitsx2sy2_so(Sx,Sy,z,p.gaussrange,p.ax);
% sxp.ztruepos=p.ztruepos;
xlim(p.ax,[p.gaussrange])
end

function fitzpar=cal_fitzpar(b,p)
 z=vertcat(b(:).z);
Sx=vertcat(b(:).sx);
Sy=vertcat(b(:).sy);
zrange=[p.gaussrange(1) p.gaussrange(2)];
fitzpar=getzfitpar(Sx,Sy,z,zrange,0,true,p.ax);
end

function out=fitsx2sy2_so(sx,sy,z,zrange,ax)

% indf=abs(z)<zrange;
indf=z>zrange(1)&z<zrange(2);

ds=sx.^2-sy.^2;
q=quantile(ds(indf),[0.05 0.95])+[-2 2];

inds=ds>q(1)&ds<q(2);
inds=indf&inds;


if sum(indf)>5
% fitpsx=fit(sx(indf).^2-sy(indf).^2,z(indf),'smoothingSpline','Normalize','on','SmoothingParam',0.95);
warning('off','curvefit:fit:iterationLimitReached');
% fitpsx=fit(sx(inds).^2-sy(inds).^2,z(inds),'poly6','Robust','LAR','Normalize','on');
fitpsx=fit(sx(inds).^2-sy(inds).^2,z(inds),'smoothingspline','Normalize','on','SmoothingParam',0.8);
indgood=fitpsx(ds)>zrange(1)&fitpsx(ds)<zrange(2);
ds2range=[min(ds(indgood)) max(ds(indgood))];
% lastwarn
% yyaxis right
% hold on
% % fitpsx=polyfit(zcorr,sx.^2-sy.^2,4);
% plot(z,sx.^2-sy.^2,'r.')
if nargin>4
    
    plot(ax,z(inds),sx(inds).^2-sy(inds).^2,'.')
    hold on
    sxsort=sort(sx.^2-sy.^2);
    zsort=feval(fitpsx,sxsort);

    plot(ax,zsort,sxsort,'k')
    plot(ax,zrange,[ds2range(1) ds2range(1)],zrange,[ds2range(2) ds2range(2)])
    ylabel(ax,'sx^2-sy^2')
     ylim(ax,[q(1) q(2)]);
    % plot(zcorr,polyval(fitpsx,zcorr),'.')

end

% xlim([-6 6])
else
    fitpsx=[];
end
out.function=fitpsx;
out.ds2range=ds2range;
end
