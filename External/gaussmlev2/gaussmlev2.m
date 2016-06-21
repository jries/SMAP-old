%gaussmlev2  MLE of single molecule positions 
%
%   [P CRLB LL t]=gaussmlev2(data,PSFSigma,iterations,fittype,Ax,Ay,Bx,By,gamma,d)
%
%   This code performs a maximum likelihood estimate of particle position,
%   emission rate (Photons/frame), and background rate (Photons/pixels/frame).
%   The found values are used to calculate the Cramer-Rao Lower Bound (CRLB)
%   for each parameter and the CRLBs are returned along with the estimated parameters.
%   The input is one or more identically sized images, each of which is
%   assumed to contain a single fluorophore.  Since the maximum likelihood
%   estimate assumes Poisson statistics, images must be converted from arbitrary
%   ADC units to photon counts.
%
%       INPUTS:
%   data:       PxPxN stack of images where N is number of PxP images
%   PSFSigma:   Microscope PSF sigma (sigma_0 for z fit, starting sigma for sigma fit)
%   iterations: Number of iterations (default=10)
%   fittype:    (default=1)
%    1: XY Position, Photons, Background
%    2: XY Position, Photons, Background, PSF Sigma
%    3: XYZ Position, Photons, Background 
%    4: XY Position, Photons, Background, PSF Sigma_X, PSF Sigma_Y
%   
%       OUTPUTS:
%   P:      Found parameter values by fittype:
%    1: [X Y Photons BG] 
%    2: [X Y Photons BG Sigma] 
%    3: [X Y Photons BG Z]  
%    4: [X Y Photons BG SigmaX SigmaY] 
%   CRLB:   Cramer-Rao lower bound calculated using P
%   LL:     Log-Likelihood calculated using Stirling's approximation      
%   t:      Execution time of fitting code. 
%
%   This codes attempts first the GPU implementation (Nvidia CUDA),
%   then the c implementation, and finally the matlab implementation
%
%   REFERENCE:
%   "Fast, single-molecule localization that
%   achieves theoretical optimal accuracy." Carlas S. Smith, Nikolai Joseph,
%   Bernd Rieger and Keith A. Lidke

function [P CRLB LL t]=gaussmlev2(data,PSFSigma,iterations,fittype,Ax,Ay,Bx,By,gamma,d)

functions = {'gaussmlev2_cuda50','gaussmlev2_cuda42','gaussmlev2_cuda40',...
    'gaussmlev2_c_thread','gaussmlev2_c','gaussmlev2_matlab'};

if nargin<3
    iterations=10;
end
if nargin<4
    fittype=1;
end
if nargin<2
    error('Minimal usage: gaussmlev2(data,PSFSigma)');
end

%convert dipimage data to single
if isa(data,'dip_image')
    data=permute(single(data),[2 1 3]);
end

%convert other matlab datatyes to single
if ~ (isa(data,'double'))
    data=single(data);
end

if ~(isa(data,'single'))
    error('Input data must be type dipimage, single, or double')
end

for ii = functions
try
tic
switch fittype
    case 1
        [P CRLB LL]=feval(ii{1},data,PSFSigma,iterations,fittype);
    case 2
        [P CRLB LL]=feval(ii{1},data,PSFSigma,iterations,fittype);
    case 3
        [P CRLB LL]=feval(ii{1},data,PSFSigma,iterations,fittype,Ax,Ay,Bx,By,gamma,d);
    case 4
        [P CRLB LL]=feval(ii{1},data,PSFSigma,iterations,fittype);    
end
t=toc;
return;
catch ME
    % just swallow the errors
end

end

% all versions must have failed 
throw ME;







