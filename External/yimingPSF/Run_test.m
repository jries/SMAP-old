%% generate experimental PSF

%07
fname = 'D:\Yiming\Data\JLM calibration\170203_olddata_161017_JB_beads_stack\1-Pos_000_000\others\image_all.tif';
saveName = 'D:\Yiming\Data\JLM calibration\170203_olddata_161017_JB_beads_stack\1-Pos_000_000\others\temp_3D-PSF_138_69_50_25_4.tif';

frz0 = 20;
pixSzIn = 138;
pixSzOut = [69];
nMol = 12;
boxSz = 25;
interpZFactor=1;
imPsfAvg= generate_psf_lut4_YL(fname,pixSzIn, pixSzOut, nMol,frz0,saveName,boxSz,interpZFactor);


%% generate Spline Coeff

%07
fname = 'D:\Yiming\Data\JLM calibration\170203_olddata_161017_JB_beads_stack\1-Pos_000_000\others\avgpsf_xy69_12mol.mat';

[np_spline coeff]=generate_psf_to_spline_YL(fname,13);


%% work with bSpline

np_psf = np_psf/max(np_psf(:));
for i = 1:size(np_psf,3)
np_psf(:,:,i) = np_psf(:,:,i)/sum(sum(np_psf(:,:,i)));
end


b3_0=bsarray(double(np_psf),'lambda',0);
b3_0_1=bsarray(double(np_psf),'lambda',0.1);%with smoothing factor 0.1


[XX,YY,ZZ]=meshgrid(1:26,1:26,1:26);
output_0 = interp3_0(b3_0,XX,YY,ZZ);
output_0_1 = interp3_0(b3_0_1,XX,YY,ZZ);

spline = Spline3D(output_0_1);
coeff = spline.coeff;





%% fit with gaussian model

Nfits=10000         %number of images to fit
bg=20;               %background fluorescence in photons/pixel/frame
Nphotons=10000;      %expected photons/frame
Npixels=13;         %linear size of fit region in pixels.
PSFsigma=1;         %PSF sigma in pixels
Ax=0;               %Aberration Terms
Bx=0;
Ay=0;
By=0;
gamma=.5 ;          %separation of x/y focal planes
d=0.5;
z=.6;               %Fixed z position
fittype=3;
%   Generate a stack of images
out = zeros(Npixels,Npixels,Nfits);
% coordsxyz = zeros(Nfits,3);
coordsx = linspace(Npixels/2-1,Npixels/2+1,10000);
coordsy = linspace(Npixels/2-1,Npixels/2+1,10000);
% coordsxy = Npixels/2 -2 +4*rand([Nfits 2]);
% coordsx = linspace(Npixels/2-1,Npixels,10000);
% coordsy = linspace(Npixels/2-1,Npixels,10000);
coordsz = linspace(-0.6,0.6,10000);
coordsxyz = [coordsx' coordsy' coordsz'];
% coordsxyz = [coordsxy coordsz'];
for i = 1:Nfits
%     coords=Npixels/2-1+rand([1 2]);
%     z = 1.2*rand(1,1)-0.6;
coords = coordsxyz(i,1:2);
z = coordsxyz(i,3);
Sx=PSFsigma*sqrt(1+((z-gamma)/d).^2+Ax*((z-gamma)/d).^3+Bx*((z-gamma)/d).^4);
Sy=PSFsigma*sqrt(1+((z+gamma)/d).^2+Ay*((z+gamma)/d).^3+By*((z+gamma)/d).^4);
out(:,:,i) = finitegausspsf(Npixels,[Sx Sy],Nphotons,bg,coords);
% coordsxyz(i,:)=[coords z];
i
end
%   Corrupt with Poisson noise
% data=single(noise(out,'poisson',1)); %requires DipImage\
data = poissrnd(out,13,13,Nfits);
[P_CS CRLB LL]=GPUmleFit_LM(single(data),single(coeff*4),200,5,0);
figure,plot(P_CS(:,end),'.')

%% test for Spline3D_v2
np_psf = double(psf);
np_psf = np_psf/max(np_psf(:));
np_psf = np_psf(25-18:25+18,25-18:25+18,:);
for i = 1:size(np_psf,3)
    np_psf(:,:,i) = np_psf(:,:,i)/sum(sum(np_psf(:,:,i)));
end
spline = Spline3D_v2(np_psf);
coeff = spline.coeff;


im =imstack(338-6:338+6,220-6:220+6,40-12:40+13);
[P] =  kernel_MLEfit_Spline_LM_SMAP_v2(im(:,:,:),single(coeff*4),13,200);

