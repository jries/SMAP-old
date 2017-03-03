function SXY=getsxyinit(p,X,Y,Z)
    SXY.Xrangeall=p.Xrange/p.cam_pixelsize_nm;
    SXY.Yrangeall=p.Yrange/p.cam_pixelsize_nm;
    SXY.Zrangeall=p.Zrange/p.cam_pixelsize_nm;
    SXY.posind=[X,Y,Z];                        
    SXY.Xrange=[p.Xrange(X), p.Xrange(X+1)]/p.cam_pixelsize_nm;
    SXY.Yrange=[p.Yrange(Y) ,p.Yrange(Y+1)]/p.cam_pixelsize_nm;
    SXY.Zrange=[p.Zrange(Z), p.Zrange(Z+1)]/p.cam_pixelsize_nm;
%     SXY.Zoffset=p.ztruepos+mean(SXY(X,Y,Z).Zrange);
    SXY.legend=[num2str(X) num2str(Y) num2str(Z)];
    SXY.Sx2_Sy2=[];
    SXY.fitzpar=[];
    SXY.splineLUT=[];
    SXY.curve=[];
    SXY.splinefit=[];
    SXY.spline=[];
    SXY.EMon=p.EMon;
end