%% This script demonstrates the use of GPUgaussMLE

clear all;

Nfits=10000  %number of images to fit
bg=1;           %background fluorescence in photons/pixel/frame
Nphotons=100;   %expected photons/frame
Npixels=7;      %linear size of fit region in pixels. 
PSFsigma=1;     %PSF sigma in pixels

%generate a stack of images
coords=Npixels/2-1+rand([Nfits 2])+[0*ones(Nfits,1) zeros(Nfits,1)];
[out] = finiteGaussPSFerf(Npixels,PSFsigma,Nphotons,bg,coords);

%corrupt with Poisson noise 
data=noise(out,'poisson',1);
iterations=5;

%fit and calculate speed
[P CRLB LL t]=mGPUgaussMLE(permute(single(data),[2 1 3]),PSFsigma,iterations,1);

CRLBx=CRLB(:,1);
X=P(:,1);

fprintf('GPUgaussMLE has performed %g fits per second.\n',Nfits/t)

%report some details
s_x_found=std(X-coords(:,1));
meanCRLBx=mean(CRLBx);

fprintf('The standard deviation of x-position error is %g \n',s_x_found)
fprintf('The mean returned CRLB based x-position uncertainty is %g \n',meanCRLBx)

%F=[CRLBx CRLBy CRLBn CRLBb]
%out1=out(:,:,0);
%[mean(N) mean(BG)]
%hist(N)