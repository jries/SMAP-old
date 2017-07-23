%% example script for fitting of 3D data
% put disclaimer

% add path to helper functions
addpath('bfmatlab')
addpath('shared')

%% make bead calibration
%run 3D calibration GUI (alternatively, you can directly call calibrate3D with proper parameters)
%and make 3D calibration
%e.g. using bead files: in /beadstacks (extracted beadstacks.zip)
% save e.g. as data/bead_3dcal.mat

calibrate3D_GUI

%% or load calibration
cal=load('data/bead_3dcal.mat'); %load bead calibration
% cal=load('data/bead_3dcal_temp.mat'); %load bead calibration



%% either simulate data or load experimental tiff file
mode =1;
if mode ==1 % simulate data
    p.offset=0;
    p.conversion=1;
    
    numlocs=1000; %number of simulated molecules. To test fitter performance on GPU, it should be >10^5
    RoiPixelsize=17; %ROI in pixels for simulation
    dz=cal.cspline.dz;  %coordinate system of spline PSF is corner based and in units pixels / planes
    z0=cal.cspline.z0; % distance and midpoint of stack in spline PSF, needed to translate into nm coordinates
    dx=floor(RoiPixelsize/2); %distance of center from edge
    ground_truth.z=linspace(-800,800,numlocs)'; %define some coordinates. Alternatively, use rand
    ground_truth.x=linspace(-0.5,0.5,numlocs)';
    ground_truth.y=sin(ground_truth.x*4*pi);
    coordinates=horzcat(ground_truth.x+dx,ground_truth.y+dx,ground_truth.z/dz+z0);
    Intensity=5000; %number of photons / localization
    background=5; %number of background photons per pixel
    imstack = simSplinePSF(RoiPixelsize,cal.cspline.coeff,Intensity,background,coordinates);
    
else %experimental data
    
    file='data/single_bead.tif'; %simulated test data, based on real bead file.
    imstackadu=readfile_ome(file); % Stack of ROIs in photons.
    p.offset=400;p.conversion=.1;
    imstack=(single(imstackadu)-p.offset)/p.conversion;% if necessary, convert ADU into photons.
    ground_truth.z=((1:size(imstack,3))'-size(imstack,3)/2+1)*10; %dz=10 nm;
end

%% fit
numlocs=size(imstack,3); %imstack needs to be in photons
imstack=single(imstack); %the fitters require the stacks in single-format;

%fit cspline, emCCD mode
tic
[Pcspline,CRLB]=mleFit_LM(imstack,5,50,single(cal.cspline.coeff));
tspline=toc;
disp(['cspline: ' num2str(numlocs/tspline) ' fits/s']);
dx=floor(size(imstack,1)/2);
cspline.x=Pcspline(:,1)-dx;cspline.y=Pcspline(:,2)-dx; %x,y in pixels 
cspline.z=(Pcspline(:,5)-cal.cspline.z0)*cal.cspline.dz;

%fit z, Gaussian model, emCCD mode
tic
[P,CRLB]=mleFit_LM(imstack,3,50,single(cal.gauss_zfit));
tgz=toc;
disp(['Gauss z: ' num2str(numlocs/tgz) ' fits/s']);
gaussz.x=P(:,1);gaussz.y=P(:,2); %x,y in pixels 
gaussz.z=P(:,5)*1000;

%fit sx,sy, Gaussian model, emCCD mode
tic
[P,CRLB]=mleFit_LM(imstack,4,50);
tgsxsy=toc;
disp(['cspline: ' num2str(numlocs/tgsxsy) ' fits/s']);

gausssxsy.x=P(:,1);gausssxsy.y=P(:,2); %x,y in pixels 
sx=P(:,5);sy=P(:,6);
gausssxsy.z=sxsy2z(sx,sy,cal.gauss_sx2_sy2); %z from sx, sy

% calculate error for central part
zrange=400;
inz=abs(ground_truth.z)<zrange;

cspline.zrel=cspline.z-ground_truth.z;
gaussz.zrel=gaussz.z-ground_truth.z;
gausssxsy.zrel=gausssxsy.z-ground_truth.z;
gausssxsy.zrel(isnan(gausssxsy.zrel))=0;
%plot fitted z vs ground-truth z
figure(101)
hold off
plot(ground_truth.z,gaussz.zrel-mean(gaussz.zrel(inz)),'r.')
hold on
plot(ground_truth.z,gausssxsy.zrel-mean(gausssxsy.zrel(inz)),'.','Color',[0.,0.8,0])
plot(ground_truth.z,cspline.zrel-mean(cspline.zrel(inz)),'b.')
gs=fit(ground_truth.z,gaussz.zrel-mean(gaussz.zrel(inz)),'smoothingspline','SmoothingParam',.00001);
gss=fit(ground_truth.z,gausssxsy.zrel-mean(gausssxsy.zrel(inz)),'smoothingspline','SmoothingParam',.001);
gz=fit(ground_truth.z,cspline.zrel-mean(cspline.zrel(inz)),'smoothingspline','SmoothingParam',.0001);

plot(ground_truth.z,ground_truth.z*0,'k','LineWidth',3)
plot(ground_truth.z,gs(ground_truth.z),'r','LineWidth',3)
plot(ground_truth.z,gss(ground_truth.z),'Color',[0.7,0.7,0],'LineWidth',3)
plot(ground_truth.z,gz(ground_truth.z),'b','LineWidth',3)

plot([-zrange -zrange],[-1000 1000],'c')
plot([zrange zrange],[-1000 1000],'c')

ss=std(cspline.zrel(inz))*10;
ylim([-ss ss])
xlabel('ground truth z')
ylabel('zfit - z (ground truth)')


cspline.dz=nanstd(ground_truth.z(inz)-cspline.z(inz));
zgauss.dz=nanstd(ground_truth.z(inz)-gaussz.z(inz));
gausssxsy.dz=nanstd(ground_truth.z(inz)-gausssxsy.z(inz));

legendtxt{4}='ground truth';
legendtxt{3}=['spline fit. error: ' num2str(cspline.dz)];
legendtxt{1}=['Gaussian z fit. error: ' num2str(zgauss.dz)];
legendtxt{2}=['Gaussian sx, sy fit. error: ' num2str(gausssxsy.dz)];
legend(legendtxt)

%calculate lateral error
cspline.dx=nanstd(ground_truth.x(inz)-cspline.x(inz));
cspline.dy=nanstd(ground_truth.y(inz)-cspline.y(inz));

dx_gaussz=nanstd(ground_truth.x(inz)-gaussz.x(inz));
dy_gaussz=nanstd(ground_truth.y(inz)-gaussz.y(inz));
gausssxsy.dx=nanstd(ground_truth.x(inz)-gausssxsy.x(inz));
gausssxsy.dy=nanstd(ground_truth.y(inz)-gausssxsy.y(inz));

%plot 3D scatter plot
numpoints=2000;
range=1:round(length(ground_truth.z)/numpoints):length(ground_truth.z);
figure(102);
hold off
scatter3(ground_truth.x,ground_truth.y,ground_truth.z);
hold on
scatter3(cspline.x,cspline.y,cspline.z);

%% or load 2D calibration
cal=load('data/bead_2dcal.mat'); %load bead calibration

ztruth = -275:50:275;
Nfits = 1000;
Nphotons = 2000;
Npixels = 17;
bg = 10;
dx = 132; %nm
dy = 132; %nm
dz=cal.cspline.dz;  %coordinate system of spline PSF is corner based and in units pixels / planes
z0=cal.cspline.z0; % distance and midpoint of stack in spline PSF, needed to translate into nm coordinates
for i = 1: length(ztruth)
    i
    coordsxy = Npixels/2 -1 +2*rand([Nfits 2]);
    coordsz = ztruth(i)/dz+z0*ones(Nfits,1);
    coordinates = [coordsxy coordsz];
    imstack = simSplinePSF(Npixels,cal.cspline.coeff,Nphotons,bg,coordinates);
    tic
    [P1 CRLB1 LL1 P2 CRLB2 LL2]=mleFit_LM(imstack,6,50,single(cal.cspline.coeff));
    tspline2D=toc;
    P = repmat(single(LL1>=LL2),1,7).*P1+repmat(single(LL1<LL2),1,7).*P2;
    CRLB = repmat(single(LL1>=LL2),1,5).*CRLB1+repmat(single(LL1<LL2),1,5).*CRLB2;
    LL = repmat(single(LL1>=LL2),1,1).*LL1+repmat(single(LL1<LL2),1,1).*LL2;
    
    z=(P(:,5)-z0).*dz;
    misassigned=sign(ztruth(i))~=sign(z);
    fractionmisassigned(i)=sum(misassigned)/length(misassigned)
    s(i)=std(z(~misassigned)-ztruth(i));
    
%     if ztruth(i)>0
%         PF=P1;
%         CRLBF=CRLB1;
%         LLF = LL1;
%     else
%         PF=P2;
%         CRLBF=CRLB2;
%         LLF = LL2;
%     end

    
    
%     misAsign(i) = sum(LL~=LLF)/length(LL);
    
    Cspline2D.x = P(:,1);
    Cspline2D.y = P(:,2);
    Cspline2D.z=(P(:,5)-z0).*dz;
    
    Cspline2D.CRLBx = CRLB(:,1);
    Cspline2D.CRLBy = CRLB(:,2);
    Cspline2D.CRLBz=CRLB(:,5);
    
    
    Cspline2D.s_x_found(i,1) = std( Cspline2D.x-coordsxy(:,1));
    Cspline2D.s_y_found(i,1) = std( Cspline2D.y-coordsxy(:,2));
    Cspline2D.s_z_found(i,1) = std( Cspline2D.z-ztruth(i));
    
    Cspline2D.meanCRLBx(i,1) = mean(sqrt(Cspline2D.CRLBx));
    Cspline2D.meanCRLBy(i,1) = mean(sqrt(Cspline2D.CRLBy));
    Cspline2D.meanCRLBz(i,1) = mean(sqrt(Cspline2D.CRLBz))*dz;
    
    
    Cspline2D.RMSEX(i,1) = sqrt(mean((Cspline2D.x-coordsxy(:,1)).^2));
    Cspline2D.RMSEY(i,1) = sqrt(mean((Cspline2D.y-coordsxy(:,2)).^2));
    Cspline2D.RMSEZ(i,1) = sqrt(mean((Cspline2D.z-ztruth(i)).^2));
    
    %Filter
    Cspline2DF.x = P(~misassigned,1);
    Cspline2DF.y = P(~misassigned,2);
    Cspline2DF.z= (P(~misassigned,5)-z0).*dz;
    
    Cspline2DF.CRLBx = CRLB(~misassigned,1);
    Cspline2DF.CRLBy = CRLB(~misassigned,2);
    Cspline2DF.CRLBz=CRLB(~misassigned,5);
    
    
    
    Cspline2DF.s_x_found(i,1) = std( Cspline2DF.x-coordsxy(~misassigned,1));
    Cspline2DF.s_y_found(i,1) = std( Cspline2DF.y-coordsxy(~misassigned,2));
    Cspline2DF.s_z_found(i,1) = std( Cspline2DF.z-ztruth(i));
    
    Cspline2DF.meanCRLBx(i,1) = mean(sqrt(Cspline2DF.CRLBx));
    Cspline2DF.meanCRLBy(i,1) = mean(sqrt(Cspline2DF.CRLBy));
    Cspline2DF.meanCRLBz(i,1) = mean(sqrt(Cspline2DF.CRLBz))*dz;
    
    
    Cspline2DF.RMSEX(i,1) = sqrt(mean((Cspline2DF.x-coordsxy(~misassigned,1)).^2));
    Cspline2DF.RMSEY(i,1) = sqrt(mean((Cspline2DF.y-coordsxy(~misassigned,2)).^2));
    Cspline2DF.RMSEZ(i,1) = sqrt(mean((Cspline2DF.z-ztruth(i)).^2));

end

figure,plot(ztruth,Cspline2D.meanCRLBx*dx,'-');
hold on
plot(ztruth,Cspline2D.meanCRLBy*dy,'-');
plot(ztruth,Cspline2D.meanCRLBz,'-');
plot(ztruth,Cspline2D.s_x_found*dx,'o');
plot(ztruth,Cspline2D.s_y_found*dy,'o');
plot(ztruth,Cspline2D.s_z_found,'o');

figure,plot(ztruth,Cspline2DF.meanCRLBx*dx,'-');
hold on
plot(ztruth,Cspline2DF.meanCRLBy*dy,'-');
plot(ztruth,Cspline2DF.meanCRLBz,'-');
plot(ztruth,Cspline2DF.s_x_found*dx,'o');
plot(ztruth,Cspline2DF.s_y_found*dy,'o');
plot(ztruth,Cspline2DF.s_z_found,'o');







%% depth-dependent calibration
% fit bead stacks in gel as shown above using the same fitter with the
% same, save fitted z-positions and frames
% z-positions should not be corrected for the refractive index mismatch. If
% they are, divide by the RI_mismatch_factor.

%load fitted coordinates of beads in gel
beads=(readtable('data/beads_in_gel_2.csv'));
% paramters for bead stack
p.dz=20;
% determine 3d corrections
zcorr=calibrate3Daberrations(beads,p);
save('data/3Daberrationcorrection.mat','zcorr') %save correction to be used later

%% correct for depth-induced errors
load('data/3Daberrationcorrection.mat') %correction
% fit your data with a 3D method of choice (example saved in
% 'data/cme_2um_site_uncorrected_sml.csv'
locs=readtable('data/cme_2um_site_uncorrected_sml.csv');

RI_mismatch_factor=.75; %refractive index mismatch correction

[img2d,xedges,yedges]=histcounts2(locs.x,locs.y,'BinWidth',[5 5]); %reconstruct superresolution image
h=fspecial('gauss',5,.6); %and blur a bit for rendering
figure(55);subplot(2,2,1);imagesc(xedges,yedges,filter2(h,img2d'));
axis equal

%%plot side view reconstruction
ypos=250; %position of slice in nm
slicewidth=60; %width of slice in nm
iny=locs.y>ypos-slicewidth/2 & locs.y<ypos+slicewidth/2; 
[img2dz,xedges,zedges]=histcounts2(locs.x(iny),locs.z(iny)*RI_mismatch_factor,'BinWidth',[10 10]);
figure(55);subplot(2,2,2);imagesc(xedges,zedges,filter2(h,img2dz'));
axis equal

% correct for aberrations and refractive index mismatch
objective_depth=2000; %position of objective above coverslip

locs.z_corr=correct_3Daberrations(zcorr,locs.z,objective_depth)*RI_mismatch_factor; %correct for aberrations
[img2dz,xedges,zedges]=histcounts2(locs.x(iny),locs.z_corr(iny),'BinWidth',[10 10]); %plot sideview
figure(55);subplot(2,2,4);imagesc(xedges,zedges,filter2(h,img2dz'));
axis equal