
%% Compiling GPUgaussMLEv2 from Matlab

%Requires:  
% 2016a or higher and parallel computing toolbox
% cuda 7.5 toolkit (restart MATLAB after install)

% Change to directory with source files- e.g. source\GPUgaussMLEv2 VS2010 noshared\GPUgaussMLEv2

%This will produce GPUgaussMLEv2.mex64 or equivalent for your platform. 

% setenv('MW_NVCC_PATH','/usr/local/CUDA/bin')

mexcuda GPUgaussMLEv2.cu mexFunction.cpp wrapper.cu



