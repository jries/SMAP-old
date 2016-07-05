function meanfitpar=meanexp(v,dq,rangev,ax)
q=myquantilefast(abs(v),[0.005,0.995],100000);
if nargin<2||isempty(dq)
    dq=min((q(2)-q(1))/1000,max(max(q(2)/1000,q(1)),q(2)/length(v)*5));
end
if nargin<3 ||isempty(rangev)
    rangev=[min(v) max(v)];
elseif length(rangev)==1
    rangev(2)=max(v);
end

n=rangev(1):dq:rangev(2);
h=histcounts(v,n);
nh=n(1:end-1)+(n(2)-n(1))/2;

d2=(n(end)-n(1));

fitfun=@(a,b,c,d,e,f,x)(a*exp(b*x)+c*exp(d*x)+e*exp(f*x));
hx=max(h);
startp=[hx/3 -1/(n(1)+d2/20) hx/3 -1/(n(1)+d2/8) hx/3 -1/( n(1)+d2)];
cf=fit(nh',h',fitfun,'Lower',[0 -abs(1/n(1)) 0 -abs(1/n(1)) 0 -abs(1/n(1))],...
    'Upper',[Inf -abs(2/n(end)) Inf -abs(2/n(end)) Inf -abs(2/n(end))],...
    'StartPoint',startp);

% meanfitpar=-(cf.c*cf.b^2+cf.a*cf.d^2)/(cf.c*cf.b^2*cf.d+cf.a*cf.b*cf.d^2);
% meanfitpar=-(cf.b*cf.d*cf.f*(cf.a/cf.b^2+cf.c/cf.d^2+cf.e/cf.f^2))/(cf.b*cf.d*cf.e+cf.b*cf.c*cf.f+cf.a*cf.d*cf.f);

meanfitpar=-(cf.a/cf.b^2+cf.c/cf.d^2+cf.e/cf.f^2)/(cf.a/cf.b+cf.c/cf.d+cf.e/cf.f);
if nargin>3&&~isempty(ax)
    axis(ax);
plot(nh,h);
hold on
plot(cf);
hold off
title(['mfit = ' num2str(meanfitpar) ', mean = ' num2str(mean(v))])
end
