function [P,CRLB,LogL]=callYimingFitter(varargin)
% 
% [P,CRLB,LogL]=GPUmleFit_LM_noInterp(varargin{:});

[P,CRLB,LogL]=GPUmleFit_LM(varargin{:},0);