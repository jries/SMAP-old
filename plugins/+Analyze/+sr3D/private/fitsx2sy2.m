function fitpsx=fitsx2sy2(sx,sy,z,zrange,ax)

% indf=abs(z)<zrange;
indf=z>zrange(1)&z<zrange(2);
if sum(indf)>5
% fitpsx=fit(sx(indf).^2-sy(indf).^2,z(indf),'smoothingSpline','Normalize','on','SmoothingParam',0.95);
fitpsx=fit(sx(indf).^2-sy(indf).^2,z(indf),'poly6','Robust','LAR','Normalize','on');
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