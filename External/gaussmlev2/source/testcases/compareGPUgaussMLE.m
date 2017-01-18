%% This script demonstrates the use of GPUgaussMLE, modified to compare
%% CUDA and MATLAB code

clear all;

Nfits=10        %number of images to fit
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
[P CRLB LL t P1 CRLB1 LL1 t1]=mGPUgaussMLEcompare(permute(single(data),[2 1 3]),PSFsigma,iterations,1);

CRLBx=CRLB(:,1);
X=P(:,1);
CRLBx1=CRLB1(:,1);
X1=P1(:,1);

fprintf('CUDA GPUgaussMLE has performed %g fits per second.\n\n',Nfits/t)

%report some details
s_x_found=std(X-coords(:,1));
meanCRLBx=mean(CRLBx);

s_x_found1=std(X1-coords(:,1));
meanCRLBx1=mean(CRLBx1);

fprintf('In CUDA The standard deviation of x-position error is %g \n',s_x_found)
fprintf('In CUDA The mean returned CRLB based x-position uncertainty is %g \n\n',meanCRLBx)

fprintf('MATLAB GPUgaussMLE has performed %g fits per second.\n\n',Nfits/t1)
fprintf('In MATLAB The standard deviation of x-position error is %g \n',s_x_found1)
fprintf('In MATLAB The mean returned CRLB based x-position uncertainty is %g \n',meanCRLBx1)

%F=[CRLBx CRLBy CRLBn CRLBb]
%out1=out(:,:,0);
%[mean(N) mean(BG)]
%hist(N)