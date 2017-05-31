%example script for fitting of 3D data
%% make bead calibration
%run 3D calibration GUI (alternatively, you can directly call calibrate3D)
%and make 3D calibration
%e.g. using bead files: xxx
%save as data/bead_3dcal.mat

calibrate3D_GUI
%% fit
file='c:xxx' %simulated test data, based on real bead file. 
imstack=readfile_ome(file); % Stack of ROIs in photons. if necessary, convert ADU into photons.

cal=load('data/bead_3dcal.mat');
coeff
zfit_par
%fit cspline, emCCD mode
[P,CRLB]=mleFit_LM(imstack,5,50,coeff);

%fit z, Gaussian model, emCCD mode
[P,CRLB]=mleFit_LM(imstack,3,50,zfit_par);

x=P(:,1);y=P(:,2); %x,y in pixels 
%correct for refractive index mismatch:
z=P(:,5)*RI_mismatch_factor;

%% depth-dependent calibration
% fit your data with a 3D method of choice (example saved in
% 'data/depth_example_sml.mat'
locs=load('data/depth_example_sml.mat');
% fit bead stacks in gel as shown above using the same fitter with the
% same, save fitted z-positions and frames
% z-positions should not be corrected for the refractive index mismatch. If
% they are, divide by the RI_mismatch_factor.
beads=load('data/beads_gel_stack_sml.mat');
% paramters
% determine 3d corrections
corr=calibrate_3Daberrations(beads,p);

% correct for aberrations and refractive index mismatch
locs.znm_corr=correct_3Daberrations(locs.znm,objective_depth)*RI_mismatch_factor;
