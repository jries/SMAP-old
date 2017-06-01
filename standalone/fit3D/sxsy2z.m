function z=sxsy2z(sx,sy,cal)
ds=sx.^2-sy.^2;
z=cal(ds);

end