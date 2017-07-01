%% example script for fitting of 3D data
% put disclaimer

% add path to helper functions
addpath('bfmatlab')
addpath('shared')
%% make bead calibration
%run 3D calibration GUI (alternatively, you can directly call calibrate3D with proper parameters)
%and make 3D calibration
%e.g. using bead files: xxx
%save as data/bead_3dcal.mat

calibrate3D_GUI
%% or load calibration
cal=load('data/bead_3dcal.mat'); %load bead calibration

%% parameters for data
p.offset=0;
p.conversion=1;
%% simulate data
numlocs=100000;
%parameteres

[imstack,ground_truth]=simulate_data(cal.spline_coeff,numlocs);
save('data/simulation.mat','imstack')
%% or load simulated data
load('data/simulation.mat')

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
[P,CRLB]=CPUmleFit_LM(imstack,5,50,single(cal.cspline_coeff));
tspline=toc;
disp(['cspline: ' num2str(numlocs/tspline) ' fits/s']);
z_cspline=(P(:,5)-cal.csplinecal.cspline.z0)*cal.csplinecal.cspline.dz;

%fit z, Gaussian model, emCCD mode
tic
[P,CRLB]=mleFit_LM(imstack,3,50,double(cal.gauss_zfit));
tgz=toc;
disp(['Gauss z: ' num2str(numlocs/tgz) ' fits/s']);
x_gaussz=P(:,1);y_gaussz=P(:,2); %x,y in pixels 
z_gaussz=P(:,5)*1000;

%fit sx,sy, Gaussian model, emCCD mode
tic
[P,CRLB]=mleFit_LM(imstack,4,50);
tgsxsy=toc;
disp(['cspline: ' num2str(numlocs/tgsxsy) ' fits/s']);

x=P(:,1);y=P(:,2); %x,y in pixels 
sx=P(:,5);sy=P(:,6);
z_gausssxsy=sxsy2z(sx,sy,cal.gauss_sx2_sy2); %z from sx, sy

% calculate error for central part


figure(101)
hold off
plot(ground_truth.z,ground_truth.z,'k')
hold on
plot(ground_truth.z,z_cspline)
plot(ground_truth.z,z_gaussz)
plot(ground_truth.z,z_gausssxsy)

legendtxt{1}='ground truth';
legendtxt{2}=['spline fit. error: ' num2str(sqrt(nanmean((ground_truth.z-z_cspline).^2)))];
legendtxt{3}=['Gaussian z fit. error: ' num2str(sqrt(nanmean((ground_truth.z-z_gaussz).^2)))];
legendtxt{4}=['Gaussian sx, sy fit. error: ' num2str(sqrt(nanmean((ground_truth.z-z_gausssxsy).^2)))];

legend(legendtxt)

%% depth-dependent calibration
% fit bead stacks in gel as shown above using the same fitter with the
% same, save fitted z-positions and frames
% z-positions should not be corrected for the refractive index mismatch. If
% they are, divide by the RI_mismatch_factor.
beads=(readtable('data/beads_in_gel_2.csv'));

% beads=load('data/beads_gel_stack_sml.mat');
% paramters
% determine 3d corrections
p.dz=20;
zcorr=calibrate3Daberrations(beads,p);
save('data/3Daberrationcorrection.mat','zcorr') %save correction to be used later

%% correct for depth-induced errors
load('data/3Daberrationcorrection.mat')
% fit your data with a 3D method of choice (example saved in
% 'data/depth_example_sml.mat'
RI_mismatch_factor=.75;
h=fspecial('gauss',5,.6);

locs=readtable('data/cme_2um_site_uncorrected_sml.csv');
[img2d,xedges,yedges]=histcounts2(locs.x,locs.y,'BinWidth',[5 5]);
figure(55);subplot(2,2,1);imagesc(xedges,yedges,filter2(h,img2d'));
axis equal
%plot slice
ypos=250;
slicewidth=60;
iny=locs.y>ypos-slicewidth/2 & locs.y<ypos+slicewidth/2;
[img2dz,xedges,zedges]=histcounts2(locs.x(iny),locs.z(iny)*RI_mismatch_factor,'BinWidth',[10 10]);
figure(55);subplot(2,2,2);imagesc(xedges,zedges,filter2(h,img2dz'));
axis equal

% correct for aberrations and refractive index mismatch
objective_depth=2000;

locs.z_corr=correct_3Daberrations(zcorr,locs.z,objective_depth)*RI_mismatch_factor;
[img2dz,xedges,zedges]=histcounts2(locs.x(iny),locs.z_corr(iny),'BinWidth',[10 10]);
figure(55);subplot(2,2,4);imagesc(xedges,zedges,filter2(h,img2dz'));
axis equal