%% This script demonstrates the use of GPUgaussMLE

clear all;

Nfits=200         %number of images to fit
bg=1;               %background fluorescence in photons/pixel/frame
Nphotons=1000;      %expected photons/frame
Npixels=7;         %linear size of fit region in pixels. 
PSFsigma=1;         %PSF sigma in pixels
Ax=0;
Bx=0;
Ay=0;
By=0;
gamma=0.5;          %separation of x/y focal planes
d=0.5;

%generate a stack of images
coords=Npixels/2-1+rand([Nfits 2])+[0*ones(Nfits,1) zeros(Nfits,1)];
z=.6;
Sx=PSFsigma*sqrt(1+((z-gamma)/d).^2+Ax*((z-gamma)/d).^3+Bx*((z-gamma)/d).^4);
Sy=PSFsigma*sqrt(1+((z+gamma)/d).^2+Ay*((z+gamma)/d).^3+By*((z+gamma)/d).^4);

[out] = finiteGaussPSFerf(Npixels,[Sx Sy],Nphotons,bg,coords);
%corrupt with Poisson noise 
data=noise(out,'poisson',1);

iterations=20;

%fit and calculate speed
tic;
[P CRLB LL]=mGPUgaussMLE(permute(single(data),[2 1 3]),PSFsigma,iterations,3,Ax,Ay,Bx,By,gamma,d);
t=toc;
CRLBx=CRLB(:,1);
X=P(:,1);
CRLBz=CRLB(:,5);
Z=P(:,5);
N=P(:,3);
BG=P(:,4);
fprintf('GPUgaussMLE has performed %g fits per second.\n',Nfits/t)

%report some details
s_x_found=std(X-coords(:,1));
meanCRLBx=mean(CRLBx);

fprintf('The standard deviation of x-position error is %g \n',s_x_found)
fprintf('The mean returned CRLB based x-position uncertainty is %g \n',meanCRLBx)

fprintf('The standard deviation of z error is %g \n',std(Z))
fprintf('The mean returned CRLB based z uncertainty is %g \n',mean(CRLBz))

%F=[CRLBx CRLBy CRLBn CRLBb]
%out1=out(:,:,0);
[mean(N) mean(BG) mean(Z)]
hist(Z)




