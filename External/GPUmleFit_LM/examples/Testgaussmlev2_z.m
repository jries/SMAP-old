%% This script demonstrates the use of GPUgaussMLE fittype 3

Nfits=1000         %number of images to fit
bg=1;               %background fluorescence in photons/pixel/frame
Nphotons=1000;      %expected photons/frame
Npixels=19;         %linear size of fit region in pixels. 
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
coords=Npixels/2-1+rand([Nfits 2])+[0*ones(Nfits,1) zeros(Nfits,1)];
Sx=PSFsigma*sqrt(1+((z-gamma)/d).^2+Ax*((z-gamma)/d).^3+Bx*((z-gamma)/d).^4);
Sy=PSFsigma*sqrt(1+((z+gamma)/d).^2+Ay*((z+gamma)/d).^3+By*((z+gamma)/d).^4);
[out] = finitegausspsf(Npixels,[Sx Sy],Nphotons,bg,coords);

%   Corrupt with Poisson noise 
data=single(noise(out,'poisson',1)); %requires DipImage
%data = poissrnd(out); %requires statistics toolbox
%   Can look at data (DipImage)
%dipshow(permute(data,[2 1 3])); %dipimage permuted 1st two dimensions

iterations=50;

%   Fit and calculate speed
%[P CRLB LL t]=gaussmlev2(data,PSFsigma,iterations,fittype,Ax,Ay,Bx,By,gamma,d);

%   Although it is recommended to call through gaussmlev2
%   specific implementations can be used and tested:

% tic;[P CRLB LL]=gaussmlev2_matlab(data,PSFsigma,iterations,fittype,Ax,Ay,Bx,By,gamma,d);t=toc
% tic;[P CRLB LL]=gaussmlev2_c(data,PSFsigma,iterations,fittype,Ax,Ay,Bx,By,gamma,d);t=toc

CRLBx=CRLB(:,1);
X=P(:,1);
CRLBz=CRLB(:,5);
Z=P(:,5);
N=P(:,3);
BG=P(:,4);
fprintf('GPUgaussMLE has performed %g fits per second.\n',Nfits/t)

%report some details
s_x_found=std(X-coords(:,1));

fprintf('The standard deviation of x-position error is %g \n',s_x_found)
fprintf('The mean returned CRLB based x-position uncertainty is %g \n',mean(sqrt(CRLBx)))

fprintf('The standard deviation of z error is %g \n',std(Z-z))
fprintf('The mean returned CRLB based z uncertainty is %g \n',mean(sqrt(CRLBz)))






