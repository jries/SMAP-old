%% This script demonstrates the use of GPUgaussMLE fittype 2

Nfits=10000     %number of images to fit
bg=1;           %background fluorescence in photons/pixel/frame
Nphotons=250;   %expected photons/frame
Npixels=7;      %linear size of fit region in pixels. 
PSFsigma=1.2;   %PSF sigma in pixels
fittype=2;

%   Generate a stack of images
coords=Npixels/2-1+rand([Nfits 2])+[0*ones(Nfits,1) zeros(Nfits,1)];
[out] = finitegausspsf(Npixels,PSFsigma,Nphotons,bg,coords);

%   Corrupt with Poisson noise 
data=single(noise(out,'poisson',1)); %requires DipImage
%data = poissrnd(out); %requires statistics toolbox
%   Can look at data (DipImage)
%dipshow(permute(data,[2 1 3])); %dipimage permuted 1st two dimensions

iterations=20;

%   Fit and calculate speed
[P CRLB LL]=GPUmleFit_LM(data,PSFsigma,20,2,0);

%   Although it is recommended to call through gaussmlev2
%   specific implementations can be used and tested:

% tic;[P CRLB LL]=gaussmlev2_matlab(data,PSFsigma,iterations,fittype);t=toc
% tic;[P CRLB LL]=gaussmlev2_c(data,[2 1 3]),PSFsigma,iterations,fittype);t=toc

% fprintf('gaussmlev2 has performed %g fits per second.\n',Nfits/t)

CRLBx=CRLB(:,1);
X=P(:,1);

CRLBs=CRLB(:,5);
S=P(:,5);

%   Report some details
s_x_found=std(X-coords(:,1));
Xmeanstd=mean(sqrt(CRLBx));

%   Report some details
s_s_found=std(S-PSFsigma);
Smeanstd=mean(sqrt(CRLBs));


fprintf('The standard deviation of x-position error is %g \n',s_x_found)
fprintf('The mean returned CRLB based x-position uncertainty is %g \n',Xmeanstd)

fprintf('The standard deviation of sigma error is %g \n',s_s_found)
fprintf('The mean returned CRLB based sigma uncertainty is %g \n',Smeanstd)

