%% This script demonstrates the use of GPUgaussMLE

clear all;

Nfits=200  %number of images to fit
bg=1;           %background fluorescence in photons/pixel/frame
Nphotons=500;   %expected photons/frame
Npixels=7;      %linear size of fit region in pixels. 
PSFsigma=1.2;     %PSF sigma in pixels

%generate a stack of images
coords=Npixels/2-1+rand([Nfits 2])+[0*ones(Nfits,1) zeros(Nfits,1)];
[out] = finiteGaussPSFerf(Npixels,PSFsigma,Nphotons,bg,coords);

%corrupt with Poisson noise 
data=noise(out,'poisson',1);
iterations=20;

%fit and calculate speed
tic;
[P CRLB LL]=mGPUgaussMLE(permute(single(data),[2 1 3]),PSFsigma,iterations,2);
t=toc;

fprintf('GPUgaussMLE has performed %g fits per second.\n',Nfits/t)
CRLBx=CRLB(:,1);
X=P(:,1);
CRLBs=CRLB(:,5);
S=P(:,5);

s_x_found=std(X-coords(:,1));
meanCRLBx=mean(CRLBx);

fprintf('The standard deviation of x-position error is %g \n',s_x_found)
fprintf('The mean returned CRLB based x-position uncertainty is %g \n',meanCRLBx)

fprintf('The standard deviation of sigma error is %g \n',std(S))
fprintf('The mean returned CRLB based sigma uncertainty is %g \n',mean(CRLBs))

%F=[CRLBx CRLBy CRLBn CRLBb]
%out1=out(:,:,0);
%[mean(N) mean(BG) mean(S)]
%hist(S)