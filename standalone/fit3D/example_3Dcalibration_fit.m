%example script for fitting of 3D data
%%
addpath('bfmatlab')
addpath('shared')
%% make bead calibration
%run 3D calibration GUI (alternatively, you can directly call calibrate3D)
%and make 3D calibration
%e.g. using bead files: xxx
%save as data/bead_3dcal.mat

% calibrate3D_GUI


%% fit
file='data/single_bead.tif'; %simulated test data, based on real bead file. 
imstack=readfile_ome(file); % Stack of ROIs in photons. 
imstack=single(imstack-400)/5;% if necessary, convert ADU into photons.

cal=load('data/bead_3dcal.mat');
numlocs=size(imstack,3);

%fit cspline, emCCD mode
tic
% [P,CRLB]=GPUmleFit_LM(imstack,5,50,single(cal.cspline_coeff));
[P,CRLB]=CPUmleFit_LM(imstack,5,50,single(cal.cspline_coeff));
% [P,CRLB]=mleFit_LM(imstack,5,50,single(cal.cspline_coeff));
tspline=toc;
disp(['cspline: ' num2str(numlocs/tspline) ' fits/s']);
zcspline=(P(:,5)-cal.csplinecal.cspline.z0)*cal.csplinecal.cspline.dz;

%fit z, Gaussian model, emCCD mode
tic
[P,CRLB]=mleFit_LM(imstack,3,50,cal.gauss_zfit);
tgz=toc;
disp(['Gauss z: ' num2str(numlocs/tgz) ' fits/s']);
x=P(:,1);y=P(:,2); %x,y in pixels 
%correct for refractive index mismatch:
zgaussz=P(:,5)*1000;

%fit sx,sy, Gaussian model, emCCD mode
tic
[P,CRLB]=mleFit_LM(imstack,4,50);
tgsxsy=toc;
disp(['cspline: ' num2str(numlocs/tgsxsy) ' fits/s']);

x=P(:,1);y=P(:,2); %x,y in pixels 
sx=P(:,5);sy=P(:,6);
zgausssxsy=sxsy2z(sx,sy,cal.gauss_sx2_sy2); %z from sx, sy

% std(z-gt) for all

figure(101)
% hold off
plot(zcspline)
hold on
plot(zgausssxsy)
plot(zgaussz)


%% depth-dependent calibration
% fit your data with a 3D method of choice (example saved in
% 'data/depth_example_sml.mat'
locs=load('data/depth_example_sml.mat');
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

% correct for aberrations and refractive index mismatch
locs.znm_corr=correct_3Daberrations(zcorr,locs.znm,objective_depth)*RI_mismatch_factor;
