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


%% parameters for data
p.offset=0;
p.conversion=1;

%% either simulate data
numlocs=1000; %number of simulated molecules. To test fitter performance on GPU, it should be >10^5
RoiPixelsize=13; %ROI in pixels for simulation
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


%% or load experimental tiff file
file='data/single_bead.tif'; %simulated test data, based on real bead file. 
imstackadu=readfile_ome(file); % Stack of ROIs in photons. 
imstack=(single(imstackadu)-p.offset)/p.conversion;% if necessary, convert ADU into photons.
ground_truth.z=((1:size(imstack,3))'-size(imstack,3)/2+1)*10; %dz=10 nm;

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
[P,CRLB]=mleFit_LM(imstack,3,50,double(cal.gauss_zfit));
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

%plot fitted z vs ground-truth z
figure(101)
hold off
plot(ground_truth.z,cspline.z)
hold on
plot(ground_truth.z,gaussz.z)
plot(ground_truth.z,gausssxsy.z)
plot([-zrange -zrange],[-1000 1000],'c')
plot([zrange zrange],[-1000 1000],'c')
plot(ground_truth.z,ground_truth.z,'k')

cspline.dz=nanstd(ground_truth.z(inz)-cspline.z(inz));
zgauss.dz=nanstd(ground_truth.z(inz)-gaussz.z(inz));
gausssxsy.dz=nanstd(ground_truth.z(inz)-gausssxsy.z(inz));

legendtxt{1}='ground truth';
legendtxt{2}=['spline fit. error: ' num2str(cspline.dz)];
legendtxt{3}=['Gaussian z fit. error: ' num2str(zgauss.dz)];
legendtxt{4}=['Gaussian sx, sy fit. error: ' num2str(gausssxsy.dz)];
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