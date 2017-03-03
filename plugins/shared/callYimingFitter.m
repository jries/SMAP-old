function [P,CRLB,LogL]=callYimingFitter(varargin)
% 
fitter=2;
switch fitter
    case 1
        [P,CRLB,LogL]=GPUmleFit_LM_noInterp(varargin{:});
    case 2
        if length(varargin)==5
            varargin{6}=1;
        end
        [P,CRLB,LogL]=GPUmleFit_LM(varargin{:});
end