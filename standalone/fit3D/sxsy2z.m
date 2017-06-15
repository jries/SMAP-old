function z=sxsy2z(sx,sy,cal)
ds=sx.^2-sy.^2; 
inrange=ds>cal.ds2range(1)*1.1&ds<cal.ds2range(2)*1.1;
z=zeros(size(ds))+NaN;
z(inrange)=cal.function(ds(inrange));
% z=cal.function(ds);
end