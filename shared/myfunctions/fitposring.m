function [x0,y0]=fitposring(x,y,R)
fh=@fixring;
xs=mean(x);ys=mean(y);
fitp=implicitfit(fh,[xs,ys],x,y,R);
x0=fitp(1);y0=fitp(2);

function err=fixring(par,x,y,R,d1,d2,d3)
err=sqrt(((x-par(1)).^2)+((y-par(2)).^2))-R;