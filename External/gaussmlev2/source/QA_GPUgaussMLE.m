% GPUgaussMLE version test
%% generate test data set
% plug in a set of standard frames later
% frames = 100;
% avgphotons = 500;
% psfsigma = 1.3;
% background = 1;
% kon = .0005;
% koff = .2;
% zoom = 10;
% [series,coords,highres,Ntot] = GenSRTestSeries('default',...
%     frames,10,10,avgphotons,20,psfsigma,background,kon,koff,zoom);

Nfits=100000  %number of images to fit
bg=10;           %background fluorescence in photons/pixel/frame
Nphotons=100;   %expected photons/frame
Npixels=7;      %linear size of fit region in pixels. 
PSFsigma=1.3;     %PSF sigma in pixels

%generate a stack of imagesc
coords=Npixels/2-1+rand([Nfits 2])+[0*ones(Nfits,1) zeros(Nfits,1)];
[out] = finiteGaussPSFerf(Npixels,PSFsigma,Nphotons,bg,0,coords);

out=noise(out,'poisson',1);
%data = permute(single(out),[2 1 3]);
data = single(out);
iterations=50;

fail_threshold = 0.00001;

%% run reference GPU gauss v1
[X Y N BG S CRLBx CRLBy CRLBn CRLBb CRLBs LL]=GPUgaussMLE(data,PSFsigma,iterations,1);
clear mex;
P_V1 = [X Y N BG];
CRLB_V1 = [CRLBx CRLBy CRLBn CRLBb];

% run repo GPUgaussMLEv2
[P_repo CRLB_repo LL]=GPUgaussMLEv2_repo(data,PSFsigma,iterations,1);
CRLB_repo=sqrt(CRLB_repo);
v1vsv2 = abs(P_V1-P_repo);
fail_v1v2 = nnz(v1vsv2 > fail_threshold);
clear mex;

%% run current 
[P_current CRLB_current LL]=GPUgaussMLEv2(data,PSFsigma,iterations,1);
CRLB_current=sqrt(CRLB_current);
v2current = abs(P_current-P_repo);
fail_current = nnz(v2current);
clear mex;

%% run default caching 64 thread limit 
[P_default CRLB_default LL]=GPUgaussMLEv2_shareddefault(data,PSFsigma,iterations,1);
CRLB_default=sqrt(CRLB_default);
v2default = abs(P_default-P_repo);
fail_v2default = nnz(v2default > fail_threshold);
clear mex;

% run favor L1 caching 64 thread limit
[P_L1 CRLB_L1 LL]=GPUgaussMLEv2_shared64(data,PSFsigma,iterations,1);
CRLB_L1=sqrt(CRLB_L1);
v2l1 = abs(P_L1-P_repo);
fail_v2l1 = nnz(v2l1 > fail_threshold);
clear mex;

%% run shared memory 192 thread limit
[P_shared192 CRLB_shared192 LL]=GPUgaussMLEv2_shared192(data,PSFsigma,iterations,1);
CRLB_shared192=sqrt(CRLB_shared192);
v2shared192 = abs(P_shared192-P_repo);
fail_v2shared192 = nnz(v2shared192 > fail_threshold);
clear mex;

%% run no shared 64 threads 
tic;
[P_noshared64 CRLB_noshared64 LL]=GPUgaussMLEv2_noshared64(data,PSFsigma,iterations,1);
CRLB_noshared64=sqrt(CRLB_noshared64);
runtime_noshared = toc;
v2noshared64 = abs(P_noshared64-P_repo);
fail_v2noshared64 = nnz(v2noshared64);
clear mex;

%% run no shared 192 threads
[P_noshared192 CRLB_noshared192 LL]=GPUgaussMLEv2_noshared192(data,PSFsigma,iterations,1);
CRLB_noshared192=sqrt(CRLB_noshared192);
v2noshared192 = abs(P_noshared192-P_repo);
fail_v2noshared192 = nnz(v2noshared192);
clear mex;

%% run simple C version
tic;
[P_c CRLB_c LL]=cGaussMLEv2(data,PSFsigma,iterations,1);
CRLB_c=sqrt(CRLB_c);
runtime_c = toc;
v2c = abs(P_c-P_repo);
fail_v2c = nnz(v2c > fail_threshold);
clear mex;

%% run C threaded version
tic;
[P_cthread CRLB_cthread LL]=cGaussMLEv2thread(data,PSFsigma,iterations,1);
CRLB_c=sqrt(CRLB_c);
runtime_cthread = toc;
v2cthread = abs(P_cthread-P_repo);
fail_v2cthread = nnz(v2cthread > fail_threshold);
clear mex;

%% check for bad localizations
max(abs(coords - P_repo(:,1:2)))
foo = gpuDevice;
foo.reset

