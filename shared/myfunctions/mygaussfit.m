function fitpar=mygaussfit(xdat,ydat,start)

% size(xdat)
% size(ydat)
% 
 oldopts=optimset('lsqcurvefit');
 oldopts=optimset(oldopts,'display','off');
newopts=optimset(oldopts,'MaxFunEvals',200,'MaxIter',200,'TolFun',1e-6,'Algorithm','levenberg-marquardt');
% newopts=optimset('MaxFunEvals',600,'MaxIter',600,'TolFun',1e-7);
% newopts=oldopts;
% fitpar=lsqcurvefit(@mygaussforfit,double(start),double(xdat),double(ydat));%,[0 0 0 0],[inf inf inf max(ydat)/100],newopts);

fitpar=lsqcurvefit(@mygaussforfit,double(start),double(xdat),double(ydat),[-inf -inf 0 -inf],[inf inf inf inf],oldopts);

% fitpar=lsqcurvefit(@mygaussforfit,double(start),double(xdat),double(ydat));