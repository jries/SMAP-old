%% This script demonstrates the use of GPUgaussMLE

clear all;

Nfits=20          %number of images to fit
bg=1;               %background fluorescence in photons/pixel/frame
Nphotons=1000;      %expected photons/frame
Npixels=19;         %linear size of fit region in pixels. 
PSFsigma=1;         %PSF sigma in pixels
Ax=0;
Bx=0;
Ay=0;
By=0;
gamma=.5 ;          %separation of x/y focal planes
d=0.5;

%generate a stack of images
coords=Npixels/2-1+rand([Nfits 2])+[0*ones(Nfits,1) zeros(Nfits,1)];
z=.5;
Sx_in=PSFsigma*sqrt(1+((z-gamma)/d).^2+Ax*((z-gamma)/d).^3+Bx*((z-gamma)/d).^4);
Sy_in=PSFsigma*sqrt(1+((z+gamma)/d).^2+Ay*((z+gamma)/d).^3+By*((z+gamma)/d).^4);
[out] = finiteGaussPSFerf(Npixels,[Sx_in Sy_in],Nphotons,bg,coords);
ndata=noise(out,'poisson');
%corrupt with Poisson noise 
data=noise(out,'poisson',1);

iterations=20;

%fit and calculate speed
tic;
[P CRLB LL]=mGPUgaussMLE(permute(single(data),[2 1 3]),PSFsigma,iterations,4);
t=toc;
CRLBx=CRLB(:,1);
CRLBy=CRLB(:,2);
X=P(:,1);
Y=P(:,2);
CRLBSx=CRLB(:,5);
CRLBSy=CRLB(:,6);
Sx=P(:,5);
Sy=P(:,6);
N=P(:,3);
BG=P(:,4);
fprintf('GPUgaussMLE has performed %g fits per second.\n',Nfits/t)

%report some details
s_x_found=std(X-coords(:,1));
s_y_found=std(Y-coords(:,2));

fprintf('The standard deviation of x-position error is %g \n',s_x_found)
fprintf('The mean returned CRLB based x-position uncertainty is %g \n',mean(CRLBx))
fprintf('The standard deviation of y-position error is %g \n',s_y_found)
fprintf('The mean returned CRLB based y-position uncertainty is %g \n',mean(CRLBy))
fprintf('The standard deviation of Sx error is %g \n',std(Sx))
fprintf('The mean returned CRLB based Sx uncertainty is %g \n',mean(CRLBSx))
fprintf('The standard deviation of Sy error is %g \n',std(Sy))
fprintf('The mean returned CRLB based Sy uncertainty is %g \n',mean(CRLBSy))

%F=[CRLBx CRLBy CRLBn CRLBb]
%out1=out(:,:,0);
[mean(N) mean(BG)  mean(Sx) mean(Sy)]
hist(X)









