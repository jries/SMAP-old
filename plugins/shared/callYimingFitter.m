function [P,CRLB,LogL]=callYimingFitter(varargin)
persistent fitter

cpufitter=@kernel_MLEfit_Spline_LM_SMAP_v2_nointerp;
gpufitter=@GPUmleFit_LM;

if isempty(fitter)
    try
        gpufitter(ones(7,'single'),1,10,1,0);
        fitter=2;
    catch
        disp('GPU didnt work, run on CPU (slow)');
        fitter=1;
    end
end

switch fitter
    case 1
        [P]=cpufitter(varargin{:});
    case 2
        varargingpu=varargin;
        if length(varargingpu)==5
            varargingpu{6}=1;
        end
        [P,CRLB,LogL]=gpufitter(varargingpu{:});
end