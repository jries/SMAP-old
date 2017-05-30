function [P,CRLB,LogL]=mleFit_LM(varargin)
% varargin:
% 1. imagestack (single)
% 2. fitmode
%   1 fix PSF
%   2 free PSF
%   3 Gauss fit z
%   4 fit PSFx, PSFy elliptical
%   5 cspline
% optional:
% 3. iterations (default=50)
% 4. paramters for fitters:
%   1. fix PSF: PSFxy sigma 
%   2. free PSF: start PSFxy
%   3. Gauss fit z: parameters: PSFx, Ax, Ay, Bx, By, gamma, d, PSFy
%   4. fit PSFx, PSFy elliptical: start PSFx, PSFy
%   5. cspline: cspline coefficients
% 5. varmap: Variance map for sCMOS. If size(varmap) ~= size(imagestack):
%       no sCMOS correction is used. Default= emCCD
% 6. silent (if 1)


%imstack, startpsf/coeff, iterations, fitmode, isemccd, hidereport

%determine of it runs on GPU, otherwise use CPU as default
persistent fitter
allfitters={@GPUmleFit_LM,@CPUmleFit_LM};
allfittersnames={'GPUmleFit_LM','CPUmleFit_LM'};
if isempty(fitter)
    for k=1:length(allfitters)
        try
            allfitters{k}(ones(7,'single'),1);
            fitter=k;
            break
        catch err
            % fitter did not work
        end
    end
    disp(['using: ' char(allfitters{fitter})]);
end


[P,CRLB,LogL]=allfitters{fitter}(varargin{:});
 clear(allfittersnames{fitter})

