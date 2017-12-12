function [P,CRLB,LogL]=mleFit_LM_global(varargin)
% varargin:
%imstack, sharedflag, iterations, spline coefficients, channelshift,
%fitmode, varmap

% 1. imagestack (single)
% 2. shared parameters

% optional:
% 3. iterations (default=50)
% 4. paramters for fitters:

%   cspline: cspline coefficients (single)
%   cspline for 2D PSF: as 5, but two fits to break asymmetry. cspline coefficients (single)
% 5. shifts dT
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

% [P,CRLB, LL] =CPUmleFit_LM_MultiChannel(fitstack,int32(sharedA),50,coeffh,single(dT),fitmode);

if nargin<7
    varmap=0; %emccd
end 
if nargin<6 || isempty(varargin{6})
    fitmode=5;
else
    fitmode=varargin{6};
end
%determine of it runs on GPU, otherwise use CPU as default
persistent fitter
allfitters={@GPUmleFit_LM_MultiChannel,@CPUmleFit_LM_MultiChannel};
allfittersnames={'GPUmleFit_LM_MultiChannel','CPUmleFit_LM_MultiChannel'};
if isempty(fitter)
    for k=1:length(allfitters)
        try
%             allfitters{k}(ones(7,7,1,2,'single'),single([1 1 1 0 0]),varargin{4},zeros(5,2));
            allfitters{k}(varargin{1:5});
            fitter=k;
            break
        catch err
            % fitter did not work
        end
    end
    disp(['using: ' char(allfitters{fitter})]);
end

if fitmode==6 %2D fit: find proper results
%     [P1,CRLB1,LogL1,P2,CRLB2,LogL2]=allfitters{fitter}(varargin{:});
%     ind1=LogL1>=LogL2;
%     ind2=LogL1<LogL2;
%     P=zeros(size(P1),'single');CRLB=zeros(size(CRLB1),'single');LogL=zeros(size(LogL1),'single');
%     P(ind1,:)=P1(ind1,:);P(ind2,:)=P2(ind2,:);
%     CRLB(ind1,:)=CRLB1(ind1,:);CRLB(ind2,:)=CRLB2(ind2,:);
%     LogL(ind1,:)=LogL1(ind1,:);LogL(ind2,:)=LogL2(ind2,:);
else
    [P,CRLB,LogL]=allfitters{fitter}(varargin{1:5});
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

