%% This script demonstrates the use of GPUgaussMLE fittype 4

Nfits=1000          %number of images to fit
bg=1;               %background fluorescence in photons/pixel/frame
Nphotons=1000;      %expected photons/frame
Npixels=19;         %linear size of fit region in pixels. 
PSFsigma=1;         %PSF sigma in pixels
fittype=4;

%   Generate a stack of images
coords=Npixels/2-1+rand([Nfits 2])+[0*ones(Nfits,1) zeros(Nfits,1)];
Sx_in=1.2;
Sy_in=1.9;
[out] = finitegausspsf(Npixels,[Sx_in Sy_in],Nphotons,bg,coords);

%   Corrupt with Poisson noise 
data=single(noise(out,'poisson',1)); %requires DipImage
%data = poissrnd(out); %requires statistics toolbox
%   Can look at data (DipImage)
%dipshow(permute(data,[2 1 3])); %dipimage permuted 1st two dimensions

iterations=20;

%   Fit and calculate speed
[P CRLB LL t]=gaussmlev2(data,PSFsigma,iterations,fittype);

%   Although it is recommended to call through gaussmlev2
%   specific implementations can be used and tested:

% tic;[P CRLB LL]=gaussmlev2_matlab(data,PSFsigma,iterations,fittype);t=toc
% tic;[P CRLB LL]=gaussmlev2_c(data,PSFsigma,iterations,fittype);t=toc

CRLBx=CRLB(:,1);
CRLBy=CRLB(:,2);
X=P(:,1);
Y=P(:,2);
CRLBSx=CRLB(:,5);
CRLBSy=CRLB(:,6);
Sx=P(:,5);
Sy=P(:,6);

fprintf('gaussmlev2 has performed %g fits per second.\n',Nfits/t)

%report some details
s_x_found=std(X-coords(:,1));
s_y_found=std(Y-coords(:,2));

fprintf('\nThe standard deviation of x-position error is %g \n',s_x_found)
fprintf('The mean returned CRLB based x-position uncertainty is %g \n\n',mean(sqrt(CRLBx)))

fprintf('The standard deviation of y-position error is %g \n',s_y_found)
fprintf('The mean returned CRLB based y-position uncertainty is %g \n\n',mean(sqrt(CRLBy)))

fprintf('The standard deviation of Sx error is %g \n',std(Sx))
fprintf('The mean returned CRLB based Sx uncertainty is %g \n\n',mean(sqrt(CRLBSx)))

fprintf('The standard deviation of Sy error is %g \n',std(Sy))
fprintf('The mean returned CRLB based Sy uncertainty is %g \n\n',mean(sqrt(CRLBSy)))











