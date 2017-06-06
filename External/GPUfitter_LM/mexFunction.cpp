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
#include <definitions.h>

#ifndef max
//! not defined in the C standard used by visual studio
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
//! not defined in the C standard used by visual studio
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

//extern void kernel_MLEFit_noshared_EMCCD_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma, const int sz, const int iterations, 
//	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits); 
//
extern void kernel_MLEFit_EMCCD_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma, const int sz, const int iterations, 
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits); 

extern void kernel_MLEFit_sigma_EMCCD_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma, const int sz, const int iterations, 
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits);

extern void kernel_MLEFit_z_EMCCD_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma_x, const float Ax, const float Ay, const float Bx, 
	const float By, const float gamma, const float d, const float PSFSigma_y, const int sz, const int iterations, 
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits);

extern void kernel_MLEFit_sigmaxy_EMCCD_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma, const int sz, const int iterations, 
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits); 

extern void kernel_splineMLEFit_z_EMCCD_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data,const float *d_coeff, const int spline_xsize, const int spline_ysize, const int spline_zsize, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits);

extern void kernel_MLEFit_sCMOS_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits, const float *d_varim); 

extern void kernel_MLEFit_sigma_sCMOS_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits,const float *d_varim); 

extern void kernel_MLEFit_z_sCMOS_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma_x, const float Ax, const float Ay, const float Bx, 
		const float By, const float gamma, const float d, const float PSFSigma_y, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits,const float *d_varim);
extern void kernel_MLEFit_sigmaxy_sCMOS_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits,const float *d_varim);



extern void kernel_splineMLEFit_z_sCMOS_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data,const float *d_coeff, const int spline_xsize, const int spline_ysize, const int spline_zsize, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits, const float *d_varim);


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
		mexErrMsgIdAndTxt("GPUmleFit_LM:cudaFail","GPUmleFit_LM(line %i): %s in %s\n",lineNumber, str, instr);
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
		mexErrMsgIdAndTxt("GPUmleFit_LM:cudaFail","%s failed with error code %i at line %d\n",str,err, lineNumber);
		exit(1); // might not stop matlab
	}
}

//*******************************************************************************************
void cudaavailable(int silent) {
	int driverVersion=0, runtimeVersion=0, deviceCount=0;
    if (cudaSuccess == cudaDriverGetVersion(&driverVersion)) {
if  (silent==0)
       mexPrintf("CUDA driver version: %d\n", driverVersion);

    } else { 
       mexErrMsgIdAndTxt("GPUmleFit_LM:nodriver","Could not query CUDA driver version\n");
    }
	if (cudaSuccess == cudaRuntimeGetVersion(&runtimeVersion)) {
if  (silent==0)
       mexPrintf("CUDA driver version: %d\n", runtimeVersion);

    } else { 
       mexErrMsgIdAndTxt("GPUmleFit_LM:noruntime","Could not query CUDA runtime version\n");
    }
	if (cudaSuccess == cudaGetDeviceCount(&deviceCount)) {
if  (silent==0)
       mexPrintf("CUDA devices detected: %d\n", deviceCount);

	} else {
       mexErrMsgIdAndTxt("GPUmleFit_LM:nodevices","Could not query CUDA device count\n", runtimeVersion);
	}
	if (deviceCount < 1) {
	   mexErrMsgIdAndTxt("GPUmleFit_LM:NoDevice","No CUDA capable devices were detected");
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
	//#define _DEBUG 1
	int blockx=0;
	int threadx=0;
	const mwSize *datasize=0;
	float PSFSigma=1, Ax=0, Ay=0, Bx=0, By=0, gamma=0.5, d=0.5, PSFSigma_y=1;
	int iterations=50;
	double * startParameters = 0;
	float *data=0, *d_data=0;
	float *d_Parameters=0,*d_CRLBs=0,*d_LogLikelihood=0;
	size_t Ndim, Nfitraw, fittype, sz, cameratype=0;//default cameratype 0 EMCCD
	mwSize noParameters;

	//sCMOS
	//float *varim=0, *d_varim=0, *gainim=0, *d_gainim=0;
	float *varim=0, *d_varim=0;

	//Spline
	int spline_xsize, spline_ysize, spline_zsize;
	float *coeff=0,*d_coeff=0;
	const mwSize *datasize_spline;

	int silent = 0;

	//if (deviceProp.kernelExecTimeoutEnabled)
	//	mexWarnMsgTxt("Warning, Kernel Execution Timeout is enabled for the GPU you are using.\nIf your fitting takes longer than the timeout it will fail.");

	//fun with timing
	LARGE_INTEGER freq;
	double togpu=0,fit=0,fromgpu=0,cleanup=0;
	LARGE_INTEGER start,stop;
	QueryPerformanceFrequency( &freq );
	QueryPerformanceCounter(&start);


	//input checks
	if (nrhs<2)
		mexErrMsgIdAndTxt("GPUmleFit_LM:WrongNoArgs","Inputs must include: GPUmleFit_LM(data,fittype,iterations,startparameter,var,silent) !\n");
	else
	{
		if (mxGetClassID(prhs[0])!=mxSINGLE_CLASS)
		mexErrMsgIdAndTxt("GPUmleFit_LM:WrongDataType","Data must be comprised of single floats!\n");

		datasize=mxGetDimensions(prhs[0]);
		Ndim=mxGetNumberOfDimensions(prhs[0]);

		if (Ndim==2)Nfitraw=1;else Nfitraw=datasize[2];

		if (datasize[0] > IMSZBIG)
			mexErrMsgIdAndTxt("GPUmleFit_LM:SubregionSize","X,Y dimension of data must be smaller than 21.\n");
		if (datasize[1]!=datasize[0])
			mexErrMsgIdAndTxt("GPUmleFit_LM:SubregionShape","Fit Box must be square");
		sz=datasize[0];

		data=(float *) mxGetData(prhs[0]);

		fittype = (int) mxGetScalar(prhs[1]);

		if (nrhs>2)
			iterations=(int) mxGetScalar(prhs[2]);

		if (nrhs>3){
			switch(fittype){
			case 1://fixed sigma
			case 2://free sigma
			case 4://sigmax and sigmay
				PSFSigma = (float)mxGetScalar(prhs[3]);
				break;
			case 3:// fit with z
				//PSFSigma = (float)mxGetScalar(prhs[3]);
				noParameters = mxGetNumberOfElements(prhs[3]);
				startParameters = mxGetPr(prhs[3]);
				if (noParameters>0)
					PSFSigma=(float) startParameters[0];
				if (noParameters>1)
					Ax = (float) startParameters[1];
				if (noParameters>2)
					Ay = (float) startParameters[2];
				if (noParameters>3)
					Bx = (float) startParameters[3];
				if (noParameters>4)
					By = (float) startParameters[4];
				if (noParameters>5)
					gamma = (float) startParameters[5];
				if (noParameters>6)
					d = (float) startParameters[6];
				if (noParameters>7)
					PSFSigma_y = (float) startParameters[7];
				else
					PSFSigma_y =PSFSigma;
				break;
			case 5://spline fit
				iterations=(int) mxGetScalar(prhs[2]);
				if (mxGetClassID(prhs[3])!=mxSINGLE_CLASS)
					mexErrMsgIdAndTxt("GPUmleFit_LM:WrongDataType","coeff must be comprised of single floats!\n");

				coeff =(float *) mxGetData(prhs[3]);
				datasize_spline=mxGetDimensions(prhs[3]);
				spline_xsize = datasize_spline[0];
				spline_ysize = datasize_spline[1];
				spline_zsize = datasize_spline[2];
				break;
			}
		}

		if (nrhs>4){
			if (mxGetNumberOfElements(prhs[4])!=mxGetNumberOfElements(prhs[0]))
				cameratype = 0;//EMCCD
			else 
			{
				if (mxGetClassID(prhs[4])!=mxSINGLE_CLASS)
					mexErrMsgIdAndTxt("GPUmleFit_LM:WrongDataType","var must be comprised of single floats!\n");
				varim=(float *) mxGetData(prhs[4]);
				cameratype = 1;//sCMOS
			}
		}

		if (nrhs>5)
			silent = (int) mxGetScalar(prhs[5]);
	}

	//query GPUs
	int deviceCount=0;
	int driverVersion=0;
	int runtimeVersion=0;
	cudaDeviceProp deviceProp; 

	cudaavailable(silent);

	cudasafe(cudaGetDeviceCount(&deviceCount),"Error detecting CUDA devices",__LINE__);
	if (deviceCount < 1)
		mexErrMsgIdAndTxt("GPUmleFit_LM:NoDevice","No CUDA capable devices were detected");

	cudasafe(cudaGetDeviceProperties(&deviceProp, 0),"Could not get properties for device 0.",__LINE__);
	cudasafe(cudaSetDevice(0),"Could not select GPU 0.",__LINE__);
	//Reset so we have a clean state before starting
	cudasafe(cudaDeviceReset(),"Error on cudaDeviceReset().",__LINE__); //release context so future cudaSetDevice calls work
	
if  (silent==0)
	mexPrintf("Using GPU %s\n",deviceProp.name);

	int BlockSize = BSZ; //for no shared scenario, we don't care!

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
	case 5:
		requiredMemory*=NV_PS;
		requiredMemory +=spline_xsize*spline_ysize*spline_zsize*64*sizeof(float);
		break;
	}

	if (cameratype ==0){//EMCCD
		requiredMemory += sz*sz*Nfitraw*sizeof(float);
	}
	else if(cameratype ==1){//sCMOS
		requiredMemory += 2*sz*sz*Nfitraw*sizeof(float);
	}
	if (requiredMemory > 0.95*deviceProp.totalGlobalMem)
		mexErrMsgIdAndTxt("GPUgaussMLEv2:NotEnoughMemory","Trying to allocation %dMB. GPU only has %dMB. Please break your fitting into multiple smaller runs.\n",
		requiredMemory/(1024*1024),availableMemory);
if  (silent==0)
{
	mexPrintf("subregion size: %d fittype: %d\n", sz, fittype);
	//mexPrintf("Number of fits: %d, Number padded fits: %d\n", Nfitraw, Nfitspad);
	mexPrintf("Number of fits: %d no padding\n", Nfitraw);
	mexPrintf("Memory required: %dMB Memory available: %dMB\n",requiredMemory/(1024*1024),availableMemory);
}

	//create device variable for data and copy to device
	cudasafe(cudaMalloc((void**)&d_data, sz*sz*Nfitraw*sizeof(float)),"Failed cudaMalloc on d_data.",__LINE__);
	cudasafe(cudaMemset(d_data, 0, sz*sz*Nfitraw*sizeof(float)),"Failed cudaMemset on d_data.",__LINE__);
	cudasafe(cudaMemcpy(d_data, data, sz*sz*Nfitraw*sizeof(float), cudaMemcpyHostToDevice),"Failed cudaMemcpy on d_data.",__LINE__);
	//sCMOS
	if(cameratype ==1){
		cudasafe(cudaMalloc((void**)&d_varim, sz*sz*Nfitraw*sizeof(float)),"Failed cudaMalloc on d_varim.",__LINE__);
		cudasafe(cudaMemset(d_varim, 0, sz*sz*Nfitraw*sizeof(float)),"Failed cudaMemset on d_varim.",__LINE__);
		cudasafe(cudaMemcpy(d_varim, varim, sz*sz*Nfitraw*sizeof(float), cudaMemcpyHostToDevice),"Failed cudaMemcpy on d_varim.",__LINE__);

		//cudasafe(cudaMalloc((void**)&d_gainim, sz*sz*Nfitraw*sizeof(float)),"Failed cudaMalloc on d_gainim.",__LINE__);
		//cudasafe(cudaMemset(d_gainim, 0, sz*sz*Nfitraw*sizeof(float)),"Failed cudaMemset on d_gainim.",__LINE__);
		//cudasafe(cudaMemcpy(d_gainim, gainim, sz*sz*Nfitraw*sizeof(float), cudaMemcpyHostToDevice),"Failed cudaMemcpy on d_gainim.",__LINE__);

	}

	if(fittype == 5){
		cudasafe(cudaMalloc((void**)&d_coeff, spline_xsize*spline_ysize*spline_zsize*64*sizeof(float)),"Failed cudaMalloc on d_coeff.",__LINE__);
		cudasafe(cudaMemset(d_coeff, 0, spline_xsize*spline_ysize*spline_zsize*64*sizeof(float)),"Failed cudaMemset on d_coeff.",__LINE__);
		cudasafe(cudaMemcpy(d_coeff, coeff, spline_xsize*spline_ysize*spline_zsize*64*sizeof(float), cudaMemcpyHostToDevice),"Failed cudaMemcpy on d_coeff.",__LINE__);
	}
	
	//create output for parameters and CRLBs
	switch(fittype){
	case 1: // (x,y,bg,I)
		cudasafe(cudaMalloc((void**)&d_Parameters,   (NV_P+1)*Nfitraw*sizeof(float)),"Failed cudaMalloc on d_Parameters",__LINE__);
		cudasafe(cudaMemset(d_Parameters, 0, (NV_P+1)*Nfitraw*sizeof(float)),"Failed cudaMemset on d_Parameters.",__LINE__);
		cudasafe(cudaMalloc((void**)&d_CRLBs,        NV_P*Nfitraw*sizeof(float)),"Failed cudaMalloc on d_CRLBs.",__LINE__);
		cudasafe(cudaMemset(d_CRLBs, 0, NV_P*Nfitraw*sizeof(float)),"Failed cudaMemset on d_CRLBs.",__LINE__);
		plhs[0]=mxCreateNumericMatrix(Nfitraw, NV_P+1, mxSINGLE_CLASS, mxREAL);
		plhs[1]=mxCreateNumericMatrix(Nfitraw, NV_P, mxSINGLE_CLASS, mxREAL);
		break;
	case 2: // (x,y,bg,I,Sigma)
	case 3: // (x,y,bg,I,z)
		cudasafe(cudaMalloc((void**)&d_Parameters,   (NV_PS+1)*Nfitraw*sizeof(float)),"Failed cudaMalloc on d_Parameters",__LINE__);
		cudasafe(cudaMemset(d_Parameters, 0, (NV_PS+1)*Nfitraw*sizeof(float)),"Failed cudaMemset on d_Parameters.",__LINE__);
		cudasafe(cudaMalloc((void**)&d_CRLBs,        NV_PS*Nfitraw*sizeof(float)),"Failed cudaMalloc on d_CRLBs.",__LINE__);
		cudasafe(cudaMemset(d_CRLBs, 0, NV_PS*Nfitraw*sizeof(float)),"Failed cudaMemset on d_CRLBs.",__LINE__);
		plhs[0]=mxCreateNumericMatrix(Nfitraw, (NV_PS+1), mxSINGLE_CLASS, mxREAL);
		plhs[1]=mxCreateNumericMatrix(Nfitraw, NV_PS, mxSINGLE_CLASS, mxREAL);
		break;      
	case 4: // (x,y,bg,I,Sx,Sy)
		cudasafe(cudaMalloc((void**)&d_Parameters,   (NV_PS2+1)*Nfitraw*sizeof(float)),"Failed cudaMalloc on d_Parameters",__LINE__);
		cudasafe(cudaMemset(d_Parameters, 0, (NV_PS2+1)*Nfitraw*sizeof(float)),"Failed cudaMemset on d_Parameters.",__LINE__);
		cudasafe(cudaMalloc((void**)&d_CRLBs,        NV_PS2*Nfitraw*sizeof(float)),"Failed cudaMalloc on d_CRLBs.",__LINE__);
		cudasafe(cudaMemset(d_CRLBs, 0, NV_PS2*Nfitraw*sizeof(float)),"Failed cudaMemset on d_CRLBs.",__LINE__);
		plhs[0]=mxCreateNumericMatrix(Nfitraw, (NV_PS2+1), mxSINGLE_CLASS, mxREAL);
		plhs[1]=mxCreateNumericMatrix(Nfitraw, NV_PS2, mxSINGLE_CLASS, mxREAL);
		break;
	case 5:
		cudasafe(cudaMalloc((void**)&d_Parameters,   (NV_PS+1)*Nfitraw*sizeof(float)),"Failed cudaMalloc on d_Parameters",__LINE__);
		cudasafe(cudaMemset(d_Parameters, 0, (NV_PS+1)*Nfitraw*sizeof(float)),"Failed cudaMemset on d_Parameters.",__LINE__);
		cudasafe(cudaMalloc((void**)&d_CRLBs,        NV_PS*Nfitraw*sizeof(float)),"Failed cudaMalloc on d_CRLBs.",__LINE__);
		cudasafe(cudaMemset(d_CRLBs, 0, NV_PS*Nfitraw*sizeof(float)),"Failed cudaMemset on d_CRLBs.",__LINE__);
		plhs[0]=mxCreateNumericMatrix(Nfitraw, (NV_PS+1), mxSINGLE_CLASS, mxREAL);
		plhs[1]=mxCreateNumericMatrix(Nfitraw, NV_PS, mxSINGLE_CLASS, mxREAL);

	//case 6:
	//	cudasafe(cudaMalloc((void**)&d_Parameters,   NV_PS*Nfitraw*sizeof(float)),"Failed cudaMalloc on d_Parameters",__LINE__);
	//	cudasafe(cudaMemset(d_Parameters, 0, NV_PS*Nfitraw*sizeof(float)),"Failed cudaMemset on d_Parameters.",__LINE__);
	//	cudasafe(cudaMalloc((void**)&d_CRLBs,        NV_PS*Nfitraw*sizeof(float)),"Failed cudaMalloc on d_CRLBs.",__LINE__);
	//	cudasafe(cudaMemset(d_CRLBs, 0, NV_PS*Nfitraw*sizeof(float)),"Failed cudaMemset on d_CRLBs.",__LINE__);
	//	plhs[0]=mxCreateNumericMatrix(Nfitraw, NV_PS, mxSINGLE_CLASS, mxREAL);
	//	plhs[1]=mxCreateNumericMatrix(Nfitraw, NV_PS, mxSINGLE_CLASS, mxREAL);
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

if  (silent==0)
	mexPrintf("threadx: %d,blockx: %d\n", threadx, blockx);


	switch(fittype) {
	case 1: //fit x,y,bg,I
		if (cameratype==0){
			kernel_MLEFit_EMCCD_wrapper(dimGrid, dimBlock, d_data, PSFSigma, (int) sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, (const int) Nfitraw);
		}
		else if (cameratype ==1){
			kernel_MLEFit_sCMOS_wrapper(dimGrid, dimBlock, d_data, PSFSigma, (int) sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, (const int) Nfitraw, d_varim);
		}
		break;

	case 2: //fit x,y,bg,I,sigma
		if (cameratype==0){
			//dim3 dimBlock(sz, sz);
   //             dim3 dimGrid(Nfitraw);
			kernel_MLEFit_sigma_EMCCD_wrapper(dimGrid, dimBlock, d_data, PSFSigma, (int) sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, (const int) Nfitraw);
		}
		else if (cameratype ==1){
			kernel_MLEFit_sigma_sCMOS_wrapper(dimGrid, dimBlock, d_data, PSFSigma, (int) sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, (const int) Nfitraw, d_varim);
		}
		break;

	case 3: //fit x,y,bg,I,z
		if (cameratype==0){
			kernel_MLEFit_z_EMCCD_wrapper(dimGrid, dimBlock, d_data, PSFSigma, Ax, Ay, Bx, By, gamma, d, PSFSigma_y, (int) sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, (const int) Nfitraw);
		}
		else if (cameratype ==1){
			kernel_MLEFit_z_sCMOS_wrapper(dimGrid, dimBlock, d_data, PSFSigma, Ax, Ay, Bx, By, gamma, d, PSFSigma_y, (int) sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, (const int) Nfitraw, d_varim);
		}
		//mexPrintf("A:%f B:%f gamma:%f d:%f \n",Ax,Bx,gamma,d);
		break;

	case 4: //fit x,y,bg,I,sigmax,sigmay
		if (cameratype==0){
			kernel_MLEFit_sigmaxy_EMCCD_wrapper(dimGrid, dimBlock, d_data, PSFSigma, (int) sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, (const int) Nfitraw);
		}
		else if (cameratype ==1){
			kernel_MLEFit_sigmaxy_sCMOS_wrapper(dimGrid, dimBlock, d_data, PSFSigma, (int) sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, (const int) Nfitraw, d_varim);
		}
		break;
	case 5:
		if (cameratype==0){
			kernel_splineMLEFit_z_EMCCD_wrapper( dimGrid,  dimBlock, d_data, d_coeff, spline_xsize, spline_ysize,  spline_zsize,  (int)sz,  iterations, 
				d_Parameters,  d_CRLBs,  d_LogLikelihood, (const int)Nfitraw);
			}
		else if (cameratype ==1){
			kernel_splineMLEFit_z_sCMOS_wrapper( dimGrid,  dimBlock, d_data, d_coeff, spline_xsize, spline_ysize,  spline_zsize,  (int)sz,  iterations, 
				d_Parameters,  d_CRLBs,  d_LogLikelihood, (const int)Nfitraw, d_varim);
			}


	}
	CUDAERROR("kernel",__LINE__);
	cudasafe(cudaDeviceSynchronize(),"sync failed",__LINE__);
	QueryPerformanceCounter(&stop);    
	fit = (double)(stop.QuadPart-start.QuadPart)/freq.QuadPart;
	QueryPerformanceCounter(&start);

	//copy to matlab output
	switch(fittype){
	case 1: // (x,y,bg,I) 
		cudasafe(cudaMemcpy((void*)mxGetData(plhs[0]), d_Parameters, (NV_P+1)*Nfitraw*sizeof(float), cudaMemcpyDeviceToHost),
			"cudaMemcpy failed for d_Parameters.",__LINE__);
		cudasafe(cudaMemcpy((void*)mxGetData(plhs[1]), d_CRLBs, NV_P*Nfitraw*sizeof(float), cudaMemcpyDeviceToHost),
			"cudaMemcpy failed for d_CRLBs.",__LINE__);
		break;
	case 2: // (x,y,bg,I,Sigma)
	case 3: // (x,y,bg,I,z)
		cudasafe(cudaMemcpy((void*)mxGetData(plhs[0]), d_Parameters, (NV_PS+1)*Nfitraw*sizeof(float), cudaMemcpyDeviceToHost),
			"cudaMemcpy failed for d_Parameters.",__LINE__);
		cudasafe(cudaMemcpy((void*)mxGetData(plhs[1]), d_CRLBs, NV_PS*Nfitraw*sizeof(float), cudaMemcpyDeviceToHost),
			"cudaMemcpy failed for d_CRLBs.",__LINE__);
		break;
	case 4: // (x,y,bg,I,Sx,Sy)
		cudasafe(cudaMemcpy((void*)mxGetData(plhs[0]), d_Parameters, (NV_PS2+1)*Nfitraw*sizeof(float), cudaMemcpyDeviceToHost),
			"cudaMemcpy failed for d_Parameters.",__LINE__);
		cudasafe(cudaMemcpy((void*)mxGetData(plhs[1]), d_CRLBs, NV_PS2*Nfitraw*sizeof(float), cudaMemcpyDeviceToHost),
			"cudaMemcpy failed for d_CRLBs.",__LINE__);
		break;
	case 5:
		cudasafe(cudaMemcpy((void*)mxGetData(plhs[0]), d_Parameters, (NV_PS+1)*Nfitraw*sizeof(float), cudaMemcpyDeviceToHost),
			"cudaMemcpy failed for d_Parameters.",__LINE__);
		cudasafe(cudaMemcpy((void*)mxGetData(plhs[1]), d_CRLBs, NV_PS*Nfitraw*sizeof(float), cudaMemcpyDeviceToHost),
			"cudaMemcpy failed for d_CRLBs.",__LINE__);
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
	cudasafe(cudaFree(d_coeff),"cudaFree failed on d_coeff.",__LINE__);
	//cudasafe(cudaFree(d_varim),"cudaFree failed on d_varim.",__LINE__);
	//cudasafe(cudaFree(d_gainim),"cudaFree failed on d_varim.",__LINE__);
	
	cudasafe(cudaDeviceReset(),"Error on cudaDeviceReset().",__LINE__); //release context so future cudaSetDevice calls work
	QueryPerformanceCounter(&stop);
if  (silent==0)
{
	mexPrintf("Memory copies to GPU %f seconds\n", togpu);
	mexPrintf("Actual fitting time %f Actual fits per second %f\n", fit,
		Nfitraw/fit);
	mexPrintf("Memory copies from GPU %f seconds\n", fromgpu);
	mexPrintf("Clean up from GPU %f seconds\n", (double)(stop.QuadPart-start.QuadPart)/freq.QuadPart);
}
	return;
}
