/*!
* \file mexFunction.cpp
* \author Keith Lidke
* \date January 10, 2010
* \brief This is a pure C file that contains only the mexFunction call and the wrappers 
* necessary to make Cuda calls.  
* [Parameters CRLBs LL]=GPUgaussMLEv2(data,PSFSigma,iterations, [fit type], [Ax], [Ay], [Bx], [By], [gamma], [d], [PSFSigma_x], [PSFSigma_y])
*/

/*! \mainpage GPUgaussMLE
*
* \section intro_sec Introduction
*
* This is the introduction.
*
* \section license_sec License
*
* \subsection step1 Step 1: Opening the box
*
* \section publication_sec Publications
* etc...
*/
#include <windows.h>
#include <tchar.h>
#include <io.h>
#pragma comment(lib, "kernel32.lib")

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>
#include <cuda_runtime.h>
#include "definitions.h"

#ifndef max
//! not defined in the C standard used by visual studio
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
//! not defined in the C standard used by visual studio
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

extern void kernel_MLEFit_noshared_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma, const int sz, const int iterations, 
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits); 

extern void kernel_MLEFit_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma, const int sz, const int iterations, 
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits); 

extern void kernel_MLEFit_sigma_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma, const int sz, const int iterations, 
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits);

extern void kernel_MLEFit_z_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma_x, const float Ax, const float Ay, const float Bx, 
	const float By, const float gamma, const float d, const float PSFSigma_y, const int sz, const int iterations, 
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits);

extern void kernel_MLEFit_sigmaxy_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma, const int sz, const int iterations, 
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits); 

void CUDAERROR(const char *instr,int lineNumber);
void cudasafe( cudaError_t err, char* str, int lineNumber);

//*******************************************************************************************
void CUDAERROR(const char *instr,int lineNumber) {
	cudaError_t errornum;
	const char *str;
	if (errornum=cudaGetLastError()) {
		//reset all cuda devices
		int deviceCount = 0;
		int ii = 0;
		cudasafe(cudaGetDeviceCount(&deviceCount),"cudaGetDeviceCount",__LINE__ ); //query number of GPUs
		for (ii = 0; ii< deviceCount;ii++) {
			cudaSetDevice(ii);
			cudaDeviceReset();
		}
		str=cudaGetErrorString(errornum);
		//mexPrintf("gpuGaussNDmle(line %i): %s in %s\n",lineNumber, str, instr);
		cudaDeviceReset();
		mexErrMsgIdAndTxt("GPUgaussMLEv2:cudaFail","gpuGaussNDmle(line %i): %s in %s\n",lineNumber, str, instr);
		exit(1); // might not stop matlab
	}
}

//*******************************************************************************************
void cudasafe( cudaError_t err, char* str, int lineNumber)
{
	if (err != cudaSuccess)
	{
		//reset all cuda devices
		int deviceCount = 0;
		int ii = 0;
		cudasafe(cudaGetDeviceCount(&deviceCount),"cudaGetDeviceCount",__LINE__ ); //query number of GPUs
		for (ii = 0; ii< deviceCount;ii++) {
			cudaSetDevice(ii);
			cudaDeviceReset();
		}
		mexErrMsgIdAndTxt("GPUgaussMLEv2:cudaFail","%s failed with error code %i at line %d\n",str,err, lineNumber);
		exit(1); // might not stop matlab
	}
}

//*******************************************************************************************
void cudaavailable() {
	int driverVersion=0, runtimeVersion=0, deviceCount=0;
    if (cudaSuccess == cudaDriverGetVersion(&driverVersion)) {
#ifdef _DEBUG
       mexPrintf("CUDA driver version: %d\n", driverVersion);
#endif 
    } else { 
       mexErrMsgIdAndTxt("GPUgaussMLEv2:nodriver","Could not query CUDA driver version\n");
    }
	if (cudaSuccess == cudaRuntimeGetVersion(&runtimeVersion)) {
#ifdef _DEBUG
       mexPrintf("CUDA driver version: %d\n", runtimeVersion);
#endif 
    } else { 
       mexErrMsgIdAndTxt("GPUgaussMLEv2:noruntime","Could not query CUDA runtime version\n");
    }
	if (cudaSuccess == cudaGetDeviceCount(&deviceCount)) {
#ifdef _DEBUG
       mexPrintf("CUDA devices detected: %d\n", deviceCount);
#endif 
	} else {
       mexErrMsgIdAndTxt("GPUgaussMLEv2:nodevices","Could not query CUDA device count\n", runtimeVersion);
	}
	if (deviceCount < 1) {
	   mexErrMsgIdAndTxt("GPUgaussMLEv2:NoDevice","No CUDA capable devices were detected");
	}
}

//*******************************************************************************************
void mexFunction(int nlhs, mxArray *plhs[],	int	nrhs, const	mxArray	*prhs[]) {
	/*!
	*  \brief Entry point in the code for Matlab.  Equivalent to main().
	*  \param nlhs number of left hand mxArrays to return
	*  \param plhs array of pointers to the output mxArrays
	*  \param nrhs number of input mxArrays
	*  \param prhs array of pointers to the input mxArrays.
	*/

	int blockx=0;
	int threadx=0;
	const mwSize *datasize=0;
	float PSFSigma=0, Ax=0, Ay=0, Bx=0, By=0, gamma=0, d=0, PSFSigma_y=0;
	int iterations=0;
	float *data=0, *d_data=0;
	float *d_Parameters=0,*d_CRLBs=0,*d_LogLikelihood=0;
	size_t Ndim, Nfitraw, fittype, sz;


	//query GPUs
	int deviceCount=0;
	int driverVersion=0;
	int runtimeVersion=0;
	cudaDeviceProp deviceProp; 

	cudaavailable();

	cudasafe(cudaGetDeviceCount(&deviceCount),"Error detecting CUDA devices",__LINE__);
	if (deviceCount < 1)
		mexErrMsgIdAndTxt("GPUgaussMLEv2:NoDevice","No CUDA capable devices were detected");

	cudasafe(cudaGetDeviceProperties(&deviceProp, 0),"Could not get properties for device 0.",__LINE__);
	cudasafe(cudaSetDevice(0),"Could not select GPU 0.",__LINE__);
	//Reset so we have a clean state before starting
	cudasafe(cudaDeviceReset(),"Error on cudaDeviceReset().",__LINE__); //release context so future cudaSetDevice calls work
	//choose larger cache
	//cudasafe(cudaDeviceSetCacheConfig(cudaFuncCachePreferL1),"Error in cudaDeviceSetCacheConfig.",__LINE__); 
	//cudasafe(cudaDeviceSetCacheConfig(cudaFuncCachePreferShared),"Error in cudaDeviceSetCacheConfig.",__LINE__); 
#ifdef _DEBUG
	mexPrintf("Using GPU %s\n",deviceProp.name);
#endif
	if (deviceProp.kernelExecTimeoutEnabled)
		mexWarnMsgTxt("Warning, Kernel Execution Timeout is enabled for the GPU you are using.\nIf your fitting takes longer than the timeout it will fail.");


	//input checks
	if (nrhs<3)
		mexErrMsgIdAndTxt("GPUgaussMLEv2:WrongNoArgs","Inputs must include: GPUgaussMLEv2(data,PSFSigma,iterations,fittype) !\n");

	if (mxGetClassID(prhs[0])!=mxSINGLE_CLASS)
		mexErrMsgIdAndTxt("GPUgaussMLEv2:WrongDataType","Data must be comprised of single floats!\n");

	datasize=mxGetDimensions(prhs[0]);
	Ndim=mxGetNumberOfDimensions(prhs[0]);

	if (Ndim==2)Nfitraw=1;else Nfitraw=datasize[2];

	if (datasize[0] > IMSZBIG)
		mexErrMsgIdAndTxt("GPUgaussMLEv2:SubregionSize","X,Y dimension of data must be smaller than 21.\n");
	if (datasize[1]!=datasize[0])
		mexErrMsgIdAndTxt("GPUgaussMLEv2:SubregionShape","Fit Box must be square");
	sz=datasize[0];
	// mexPrintf("c sizeX: %d Ndims: %d\n",sz,Ndim);

	//put as many images as fit into a block
	//float szsqaure = sz*sz;
	//int BlockSize = (int) floor((float)15000/4/sz/sz);
	//int BlockSize = (int) floor((float)MEM/(sz*sz));
	//int BlockSize = (int) floor((float)(deviceProp.sharedMemPerBlock-1024)/(sizeof(float)*sz*sz)); // assumes allocating all of shared
	//mexPrintf("BlockSize %d = (int) floor((float)(%d-1024)/sizeof(%d)*%d*%d\n",BlockSize,deviceProp.sharedMemPerBlock,sizeof(float),sz,sz);
	//BlockSize = max(4, BlockSize);
	//BlockSize = min(BSZ, BlockSize);
	int BlockSize = BSZ; //for no shared scenario, we don't care!


	//int Nfitspad= BlockSize*(int) ceil( (float) Nfitraw/BlockSize);
	//const int Nfitspad=(int) Nfitraw; //kernel auto aborts when thread > Nfitspad?

	//fun with timing
	LARGE_INTEGER freq;
	double togpu=0,fit=0,fromgpu=0,cleanup=0;
	LARGE_INTEGER start,stop;
	QueryPerformanceFrequency( &freq );
	QueryPerformanceCounter(&start);

	//get variables
	data=(float *) mxGetData(prhs[0]);
	PSFSigma=(float)mxGetScalar(prhs[1]); //matlab-dip_image convention
	iterations=(int) mxGetScalar(prhs[2]);
	fittype = (int) mxGetScalar(prhs[3]);

	if (nrhs>4)
		Ax = (float) mxGetScalar(prhs[4]);
	if (nrhs>5)
		Ay = (float) mxGetScalar(prhs[5]);
	if (nrhs>6)
		Bx = (float) mxGetScalar(prhs[6]);
	if (nrhs>7)
		By = (float) mxGetScalar(prhs[7]);
	if (nrhs>8)
		gamma = (float) mxGetScalar(prhs[8]);
	if (nrhs>9)
		d = (float) mxGetScalar(prhs[9]);
	if (nrhs>10)
		PSFSigma_y = (float) mxGetScalar(prhs[10]);
	else
		PSFSigma_y =PSFSigma;

	//check if we are doing something silly and allocating more memory than the card has
	const size_t availableMemory = deviceProp.totalGlobalMem/(1024*1024);
	size_t requiredMemory = 3*Nfitraw;
	switch (fittype) {
	case 1:
		requiredMemory*=NV_P;
		break;
	case 2:
	case 3:
		requiredMemory*=NV_PS;
		break;
	case 4:
		requiredMemory*=NV_PS2;
		break;
	}
	requiredMemory += sz*sz*Nfitraw*sizeof(float);
	if (requiredMemory > 0.95*deviceProp.totalGlobalMem)
		mexErrMsgIdAndTxt("GPUgaussMLEv2:NotEnoughMemory","Trying to allocation %dMB. GPU only has %dMB. Please break your fitting into multiple smaller runs.\n",
		requiredMemory/(1024*1024),availableMemory);
#ifdef _DEBUG
	mexPrintf("subregion size: %d fittype: %d\n", sz, fittype);
	//mexPrintf("Number of fits: %d, Number padded fits: %d\n", Nfitraw, Nfitspad);
	mexPrintf("Number of fits: %d no padding\n", Nfitraw);
	mexPrintf("Memory required: %dMB Memory available: %dMB\n",requiredMemory/(1024*1024),availableMemory);
#endif

	//create device variable for data and copy to device
	cudasafe(cudaMalloc((void**)&d_data, sz*sz*Nfitraw*sizeof(float)),"Failed cudaMalloc on d_data.",__LINE__);
	cudasafe(cudaMemset(d_data, 0, sz*sz*Nfitraw*sizeof(float)),"Failed cudaMemset on d_data.",__LINE__);
	cudasafe(cudaMemcpy(d_data, data, sz*sz*Nfitraw*sizeof(float), cudaMemcpyHostToDevice),"Failed cudaMemcpy on d_data.",__LINE__);

	//create output for parameters and CRLBs
	switch(fittype){
	case 1: // (x,y,bg,I)
		cudasafe(cudaMalloc((void**)&d_Parameters,   NV_P*Nfitraw*sizeof(float)),"Failed cudaMalloc on d_Parameters",__LINE__);
		cudasafe(cudaMemset(d_Parameters, 0, NV_P*Nfitraw*sizeof(float)),"Failed cudaMemset on d_Parameters.",__LINE__);
		cudasafe(cudaMalloc((void**)&d_CRLBs,        NV_P*Nfitraw*sizeof(float)),"Failed cudaMalloc on d_CRLBs.",__LINE__);
		cudasafe(cudaMemset(d_CRLBs, 0, NV_P*Nfitraw*sizeof(float)),"Failed cudaMemset on d_CRLBs.",__LINE__);
		plhs[0]=mxCreateNumericMatrix(Nfitraw, NV_P, mxSINGLE_CLASS, mxREAL);
		plhs[1]=mxCreateNumericMatrix(Nfitraw, NV_P, mxSINGLE_CLASS, mxREAL);
		break;
	case 2: // (x,y,bg,I,Sigma)
	case 3: // (x,y,bg,I,z)
		cudasafe(cudaMalloc((void**)&d_Parameters,   NV_PS*Nfitraw*sizeof(float)),"Failed cudaMalloc on d_Parameters",__LINE__);
		cudasafe(cudaMemset(d_Parameters, 0, NV_PS*Nfitraw*sizeof(float)),"Failed cudaMemset on d_Parameters.",__LINE__);
		cudasafe(cudaMalloc((void**)&d_CRLBs,        NV_PS*Nfitraw*sizeof(float)),"Failed cudaMalloc on d_CRLBs.",__LINE__);
		cudasafe(cudaMemset(d_CRLBs, 0, NV_PS*Nfitraw*sizeof(float)),"Failed cudaMemset on d_CRLBs.",__LINE__);
		plhs[0]=mxCreateNumericMatrix(Nfitraw, NV_PS, mxSINGLE_CLASS, mxREAL);
		plhs[1]=mxCreateNumericMatrix(Nfitraw, NV_PS, mxSINGLE_CLASS, mxREAL);
		break;      
	case 4: // (x,y,bg,I,Sx,Sy)
		cudasafe(cudaMalloc((void**)&d_Parameters,   NV_PS2*Nfitraw*sizeof(float)),"Failed cudaMalloc on d_Parameters",__LINE__);
		cudasafe(cudaMemset(d_Parameters, 0, NV_PS2*Nfitraw*sizeof(float)),"Failed cudaMemset on d_Parameters.",__LINE__);
		cudasafe(cudaMalloc((void**)&d_CRLBs,        NV_PS2*Nfitraw*sizeof(float)),"Failed cudaMalloc on d_CRLBs.",__LINE__);
		cudasafe(cudaMemset(d_CRLBs, 0, NV_PS2*Nfitraw*sizeof(float)),"Failed cudaMemset on d_CRLBs.",__LINE__);
		plhs[0]=mxCreateNumericMatrix(Nfitraw, NV_PS2, mxSINGLE_CLASS, mxREAL);
		plhs[1]=mxCreateNumericMatrix(Nfitraw, NV_PS2, mxSINGLE_CLASS, mxREAL);
		break;
	}
	cudasafe(cudaMalloc((void**)&d_LogLikelihood, Nfitraw*sizeof(float)),"Failed cudaMalloc on d_LogLikelihood.",__LINE__);
	plhs[2]=mxCreateNumericMatrix(Nfitraw, 1, mxSINGLE_CLASS, mxREAL);

	//check allocations
	if (plhs[0] == NULL || plhs[1] == NULL || plhs[2] == NULL)
		mexErrMsgIdAndTxt("GPUgaussMLEv2:NotEnoughMemory","Could not allocate memory for results.\n");

	QueryPerformanceCounter(&stop);
	togpu = (double) (stop.QuadPart-start.QuadPart)/freq.QuadPart;
	QueryPerformanceCounter(&start);

	//setup kernel
	blockx = (int) ceil( (float)Nfitraw/(float)BlockSize);
	threadx= BlockSize;

	dim3 dimBlock(threadx);
	dim3 dimGrid(blockx);

#ifdef _DEBUG
	mexPrintf("threadx: %d,blockx: %d\n", threadx, blockx);
#endif

	switch(fittype) {
	case 1: //fit x,y,bg,I
		kernel_MLEFit_noshared_wrapper(dimGrid, dimBlock, d_data, PSFSigma, (int) sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, (const int) Nfitraw);
		//kernel_MLEFit_wrapper(dimGrid, dimBlock, d_data, PSFSigma, (int) sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfitraw);
		break;

	case 2: //fit x,y,bg,I,sigma
		kernel_MLEFit_sigma_wrapper(dimGrid, dimBlock, d_data, PSFSigma, (int) sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, (const int) Nfitraw);
		break;

	case 3: //fit x,y,bg,I,z
		kernel_MLEFit_z_wrapper(dimGrid, dimBlock, d_data, PSFSigma, Ax, Ay, Bx, By, gamma, d, PSFSigma_y, (int) sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, (const int) Nfitraw);
		//mexPrintf("A:%f B:%f gamma:%f d:%f \n",Ax,Bx,gamma,d);
		break;

	case 4: //fit x,y,bg,I,sigmax,sigmay
		kernel_MLEFit_sigmaxy_wrapper(dimGrid, dimBlock, d_data, PSFSigma, (int) sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, (const int) Nfitraw);

	}
	CUDAERROR("kernel",__LINE__);
	cudasafe(cudaDeviceSynchronize(),"sync failed",__LINE__);
	QueryPerformanceCounter(&stop);    
	fit = (double)(stop.QuadPart-start.QuadPart)/freq.QuadPart;
	QueryPerformanceCounter(&start);

	//copy to matlab output
	switch(fittype){
	case 1: // (x,y,bg,I) 
		cudasafe(cudaMemcpy((void*)mxGetData(plhs[0]), d_Parameters, NV_P*Nfitraw*sizeof(float), cudaMemcpyDeviceToHost),
			"cudaMemcpy failed for d_Parameters.",__LINE__);
		cudasafe(cudaMemcpy((void*)mxGetData(plhs[1]), d_CRLBs, NV_P*Nfitraw*sizeof(float), cudaMemcpyDeviceToHost),
			"cudaMemcpy failed for d_CRLBs.",__LINE__);
		break;
	case 2: // (x,y,bg,I,Sigma)
	case 3: // (x,y,bg,I,z)
		cudasafe(cudaMemcpy((void*)mxGetData(plhs[0]), d_Parameters, NV_PS*Nfitraw*sizeof(float), cudaMemcpyDeviceToHost),
			"cudaMemcpy failed for d_Parameters.",__LINE__);
		cudasafe(cudaMemcpy((void*)mxGetData(plhs[1]), d_CRLBs, NV_PS*Nfitraw*sizeof(float), cudaMemcpyDeviceToHost),
			"cudaMemcpy failed for d_CRLBs.",__LINE__);
		break;
	case 4: // (x,y,bg,I,Sx,Sy)
		cudasafe(cudaMemcpy((void*)mxGetData(plhs[0]), d_Parameters, NV_PS2*Nfitraw*sizeof(float), cudaMemcpyDeviceToHost),
			"cudaMemcpy failed for d_Parameters.",__LINE__);
		cudasafe(cudaMemcpy((void*)mxGetData(plhs[1]), d_CRLBs, NV_PS2*Nfitraw*sizeof(float), cudaMemcpyDeviceToHost),
			"cudaMemcpy failed for d_CRLBs.",__LINE__);
		break;
	}
	cudasafe(cudaMemcpy((void*)mxGetData(plhs[2]), d_LogLikelihood, Nfitraw*sizeof(float), cudaMemcpyDeviceToHost),
		"cudaMemcpy failed for d_LogLikelihood.",__LINE__);

	QueryPerformanceCounter(&stop);
	fromgpu = (double)(stop.QuadPart-start.QuadPart)/freq.QuadPart;
	QueryPerformanceCounter(&start);

	//cleanup
	cudasafe(cudaFree(d_data),"cudaFree failed on d_data.",__LINE__);
	cudasafe(cudaFree(d_Parameters),"cudaFree failed on d_Parameters.",__LINE__);
	cudasafe(cudaFree(d_CRLBs),"cudaFree failed on d_CRLBs.",__LINE__);
	cudasafe(cudaFree(d_LogLikelihood),"cudaFree failed on d_LogLikelihood.",__LINE__);

	cudasafe(cudaDeviceReset(),"Error on cudaDeviceReset().",__LINE__); //release context so future cudaSetDevice calls work
	QueryPerformanceCounter(&stop);
#ifdef _DEBUG
	mexPrintf("Memory copies to GPU %f seconds\n", togpu);
	mexPrintf("Actual fitting time %f Actual fits per second %f\n", fit,
		Nfitraw/fit);
	mexPrintf("Memory copies from GPU %f seconds\n", fromgpu);
	mexPrintf("Clean up from GPU %f seconds\n", (double)(stop.QuadPart-start.QuadPart)/freq.QuadPart);
#endif
	return;
}
