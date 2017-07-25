function [P,CRLB,LogL,P1,CRLB1,LogL1,P2,CRLB2,LogL2]=mleFit_LM(varargin)
% varargin:
%imstack, startpsf/coeff, iterations, fitmode, isemccd, hidereport
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
%       (single)
%   4. fit PSFx, PSFy elliptical: start PSFx, PSFy
%   5. cspline: cspline coefficients (single)
%   6. cspline for 2D PSF: as 5, but two fits to break asymmetry. cspline coefficients (single)
% 5. varmap: Variance map for sCMOS. If size(varmap) ~= size(imagestack):
%       no sCMOS correction is used. Default= emCCD
% 6. silent (suppress output if 1)

%Output:
%P
%1. X, Y, Photons, Background, Iterations
%2. X, Y, Photons, Background, PSFxy, Iterations
%3. X, Y, Photons, Background, Z, Iterations
%4. X, Y, Photons, Background, PSFx, PSFy, Iterations
%5. X, Y, Photons, Background, Z, Iterations
%6. X, Y, Photons, Background, Z, Iterations
%CRLB: cramer-rao lower bounds, as in P
%LogL: log-likelihood.

%Only for fitmode 6: P1 etc: results with z-startparameter<0, P2 etc:
%results with z-startparameter>0



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

if varargin{2}==6 %2D fit: find proper results
    [P1,CRLB1,LogL1,P2,CRLB2,LogL2]=allfitters{fitter}(varargin{:});
    ind1=LogL1>=LogL2;
    ind2=LogL1<LogL2;
    P=zeros(size(P1),'single');CRLB=zeros(size(CRLB1),'single');LogL=zeros(size(LogL1),'single');
    P(ind1,:)=P1(ind1,:);P(ind2,:)=P2(ind2,:);
    CRLB(ind1,:)=CRLB1(ind1,:);CRLB(ind2,:)=CRLB2(ind2,:);
    LogL(ind1,:)=LogL1(ind1,:);LogL(ind2,:)=LogL2(ind2,:);
else
    [P,CRLB,LogL]=allfitters{fitter}(varargin{:});
    P1=[];CRLB1=[];LogL1=[];P2=[];CRLB2=[];LogL2=[];
end
%%
% 
% <latex>
% \begin{tabular}{|c|c|} \hline
% $n$ & $n!$ \\ \hline
% 1 & 1 \\
% 2 & 2 \\
% 3 & 6 \\ \hline
% \end{tabular}
% </latex>
% 
 clear(allfittersnames{fitter})

