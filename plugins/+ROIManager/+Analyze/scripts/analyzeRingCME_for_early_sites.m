%load sites in ROImanager, evaluate with CME2DRing and intensityTiff

global se;

range=1:300;
sites=se.sites(range);


GFP1=getFieldAsVector(sites,'evaluation.intensityTiff.Amplitude1');
asymmetry1=getFieldAsVectorInd(sites,'evaluation.CME2DRing.imfit.asymmetry',1);
asymmetry2=getFieldAsVectorInd(sites,'evaluation.CME2DRing.imfit.asymmetry',2);
area=getFieldAsVector(sites,'evaluation.CME2DRing.imfit.areastat.Area');
ecc=getFieldAsVector(sites,'evaluation.CME2DRing.imfit.areastat.Eccentricity');
fractionCircle=getFieldAsVector(sites,'evaluation.CME2DRing.imfit.areastat.fractionCircle');
theta_varnorm=getFieldAsVector(sites,'evaluation.CME2DRing.imfit.profiles1.theta_varnorm');
xcorramp=getFieldAsVector(sites,'evaluation.CME2DRing.imfit.profiles1.xcorramp');
r1=getFieldAsVector(sites,'evaluation.CME2DRing.imfit.r1');
dr_ro1=getFieldAsVector(sites,'evaluation.CME2DRing.imfit.dr_ro1');
images=getFieldAsVector(sites,'image','image');



figure(9312);
clf;

mode='median';
numberofbins=10;
GFPmaxcutoff=0.5; % factor*max(GFP) to exlude brightest GFP values, which are usually sparse and artifical (overlap)
binwidth=[];
ringfraction=[];
binmidpoints=[];
ringthreshold=1.5;
gfprange=linspace(1,GFPmaxcutoff*max(GFP1),numberofbins);
for i=1:(numberofbins - 1)
    binwidth(i)=sum(GFP1>gfprange(i) & GFP1<=gfprange(i+1));
    ringfraction(i)=sum(dr_ro1(GFP1>gfprange(i) & GFP1<=gfprange(i+1))<ringthreshold)/binwidth(i);
    binmidpoints(i)=mean(gfprange(i:i+1));   
end

    

figure(9312);
subplot(3,3,1)
hold off
plot(GFP1,asymmetry1./asymmetry2,'.');
hold on
plot(gfprange,bindata(GFP1,asymmetry1./asymmetry2,gfprange,mode))
title('GFP vs asymmetry1/2');

subplot(3,3,2)
hold off
plot(GFP1,ecc,'.');
hold on
plot(gfprange,bindata(GFP1,ecc,gfprange,mode))
title('GFP vs ecc');

subplot(3,3,3)
hold off
plot(GFP1,fractionCircle,'.');
hold on
plot(gfprange,bindata(GFP1,fractionCircle,gfprange,mode))
title('GFP vs fractionCircle');

subplot(3,3,4)
hold off
plot(GFP1,theta_varnorm,'.');
hold on
plot(gfprange,bindata(GFP1,theta_varnorm,gfprange,mode))
title('GFP vs theta_varnorm');

subplot(3,3,5)
hold off
plot(GFP1,xcorramp,'.');
hold on
plot(gfprange,bindata(GFP1,xcorramp,gfprange,mode))
title('GFP vs xcorramp');

subplot(3,3,6)
hold off
plot(GFP1,area,'.');
hold on
plot(gfprange,bindata(GFP1,area,gfprange,mode))
title('GFP vs area');

subplot(3,3,7)
hold off
plot(GFP1,r1,'.');
hold on
plot(gfprange,bindata(GFP1,r1,gfprange,mode))
title('GFP vs outer radius');

subplot(3,3,8)
plot(gfprange,bindata(GFP1,dr_ro1,gfprange,mode))
title('GFP vs dr/ro_1');

subplot(3,3,9)
plot(binmidpoints,ringfraction);
title('GFP vs ring fraction');

