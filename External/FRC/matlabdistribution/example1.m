%Example1
% -- needs the matlab toolbox dipimage, free download at www.diplib.org
% run each section separtely. This file shows a lot of different
% possibilities to compute FIRE from positions via the function
% "postoresolution"
% -- you need to have the fullpaht to the FIREfunctions directory on the MATLAB path  

clear all
close all

%% loading the data
% localizations of single emitters in the format:
%  x, y, t [pixels, pixels, frames]

fprintf('\n -- loading the data --\n')
coords = dlmread(['..' filesep 'ExampleData' filesep 'example_Fig2a.dat'],',');
pixelsize = 16e3/150;                                                           %in nanometers

%% make an overview image
superzoom = 10;
szx = superzoom * 256;
szy = superzoom * 256;

im = binlocalizations(coords, szx, szy, superzoom);
h=dipshow(im);
dipmapping(h,[0 5],'colormap',hot)


%% compute FIRE
fprintf('\n -- computing FIRE --\n')
[fire_value, ~, fireH, fireL] = postoresolution(coords, szx, superzoom); % in super-resolution pixels
fprintf('FIRE value %2.1f +- %2.2f [px]\n', fire_value, (fireL-fireH)/2);
fprintf('FIRE value %2.1f +- %2.2f [nm]\n', fire_value*pixelsize/superzoom, (fireL-fireH)/2*pixelsize/superzoom);


%% compute FIRE as a function of frame time (takes ~1-2min.)
fprintf('\n -- computing FIRE as a function of time for in 25 steps--\n')
tfrac = 25;
[~,~,~,~,fireT] = postoresolution(coords, szx, superzoom, 500, tfrac); % in super-resolution pixels
figure;
plot(linspace(0,1,tfrac),fireT*pixelsize/superzoom,'x-')
xlabel('time fraction')
ylabel('FIRE (nm)')
title('FIRE as a function of total number of frames')


%% compute FRC curve
fprintf('\n -- computing FRC curve--\n')
[~,frc_curve] = postoresolution(coords, szx, superzoom); 
figure;
qmax = 0.5/(pixelsize/superzoom);
plot(linspace(0,qmax*sqrt(2), length(frc_curve)), frc_curve,'-')
xlim([0,qmax])
hold on
plot([0 qmax],[1/7 1/7],'r-');
plot([0 qmax],[0 0],'k--'); hold off
xlabel('spatial frequency (nm^{-1})')
ylabel('FRC')
title('Fourier Ring Correlation curve')


%% compute FIRE (averaged over 20 runs)
fprintf('\n -- computing FIRE averaged over 20 random block splits --\n')
[fire_value, ~, fireH, fireL] = postoresolution(coords, szx, superzoom, 500,[], 20); % in super-resolution pixels
fprintf('FIRE value %2.1f +- %2.2f [px]\n', fire_value, (fireL-fireH)/2);
fprintf('FIRE value %2.1f +- %2.2f [nm]\n', fire_value*pixelsize/superzoom, (fireL-fireH)/2*pixelsize/superzoom);

%% compute FIRE as a function of frame time (averaged over 20 runs) - takes a while 10-20min.
% uncomment to run
% fprintf('\n -- computing FIRE versus time; each time point averaged over 20 random block splits--\n')
% [~,~,~,~,fireT] = postoresolution(coords, szx, superzoom, 500, tfrac, 20); % in super-resolution pixels
% figure;
% plot(linspace(0,1,tfrac),fireT*pixelsize/superzoom,'x-')
% xlabel('time fraction')
% ylabel('FIRE (nm)')
% title('FIRE as a function of total number of frames')
