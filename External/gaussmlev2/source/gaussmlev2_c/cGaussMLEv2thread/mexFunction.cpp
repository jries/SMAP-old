/*!
* \file mexFunction.cpp
* \author Keith Lidke
* \date January 10, 2010
* \brief This is a pure C file that contains only the mexFunction call and the wrappers 
* necessary to make Cuda calls.  
* [Parameters CRLBs LL]=CgaussMLEv2(data,PSFSigma,iterations, [fit type], [Ax], [Ay], [Bx], [By], [gamma], [d], [PSFSigma_x], [PSFSigma_y])
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
#pragma comment(lib, "kernel32.lib")

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>
#include "definitions.h"
#include "GPUGaussMLEv2.h"
#include "FitInfo.h"
#include <mutex>
#include <thread>

//*******************************************************************************************
void worker(const size_t id, FitInfo &info) {
	
	int subregion=0; 

	info.worker_mutex->lock();
	while (info.current < info.Nfitraw) {
		subregion = (int) info.current++;
		info.worker_mutex->unlock();
		//mexPrintf("thread %d processing %d of %d\n",id,info.current,info.Nfitraw);

	switch(info.fittype) {
		case 1: //fit x,y,bg,I
			kernel_MLEFit(subregion, info.data, info.PSFSigma, (int) info.sz, info.iterations, info.Parameters, info.CRLBs, 
			info.LogLikelihood, (const int) info.Nfitraw);
		break;

		case 2: //fit x,y,bg,I,sigma
			kernel_MLEFit_sigma(subregion, info.data, info.PSFSigma, (int) info.sz, info.iterations, info.Parameters, 
			info.CRLBs, info.LogLikelihood, (const int) info.Nfitraw);
		break;

		case 3: //fit x,y,bg,I,z
			kernel_MLEFit_z(subregion, info.data, info.PSFSigma, info.Ax, info.Ay, info.Bx, info.By, info.gamma, info.d, info.PSFSigma_y, 
			(int) info.sz, info.iterations, info.Parameters, info.CRLBs, info.LogLikelihood, (const int) info.Nfitraw);
		//mexPrintf("A:%f B:%f gamma:%f d:%f \n",Ax,Bx,gamma,d);
		break;

		case 4: //fit x,y,bg,I,sigmax,sigmay
			kernel_MLEFit_sigmaxy(subregion, info.data, info.PSFSigma, (int) info.sz, info.iterations, info.Parameters, info.CRLBs, 
			info.LogLikelihood, (const int) info.Nfitraw);
		break;
		}

		info.worker_mutex->lock();
	}
	info.worker_mutex->unlock();
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

	//int blockx=0;
	//int threadx=0;
	const mwSize *datasize=0;
	//float PSFSigma=0, Ax=0, Ay=0, Bx=0, By=0, gamma=0, d=0, PSFSigma_y=0;
	//int iterations=0;
	//float *data=0;
	//float *Parameters=0,*CRLBs=0,*LogLikelihood=0;
	size_t Ndim;
	//, Nfitraw, fittype, sz;

	FitInfo info;

	//input checks
	if (nrhs<3)
		mexErrMsgIdAndTxt("cGaussMLEv2:WrongNoArgs","Inputs must include: cGaussMLEv2(data,PSFSigma,iterations,fittype) !\n");

	if (mxGetClassID(prhs[0])!=mxSINGLE_CLASS)
		mexErrMsgIdAndTxt("cGaussMLEv2:WrongDataType","Data must be comprised of single floats!\n");

	datasize=mxGetDimensions(prhs[0]);
	Ndim=mxGetNumberOfDimensions(prhs[0]);

	if (Ndim==2)info.Nfitraw=1;else info.Nfitraw=datasize[2];

	if (datasize[0] > IMSZBIG)
		mexErrMsgIdAndTxt("cGaussMLEv2:SubregionSize","X,Y dimension of data must be smaller than 21.\n");
	if (datasize[1]!=datasize[0])
		mexErrMsgIdAndTxt("cGaussMLEv2:SubregionShape","Fit Box must be square");
	info.sz=datasize[0];
	// mexPrintf("c sizeX: %d Ndims: %d\n",sz,Ndim);

	//query CPU cores

	//put as many images as fit into a block
	//float szsqaure = sz*sz;
	//int BlockSize = (int) floor((float)15000/4/sz/sz);
	//int BlockSize = (int) floor((float)MEM/(sz*sz));
	//int BlockSize = (int) floor((float)(deviceProp.sharedMemPerBlock-1024)/(sizeof(float)*sz*sz)); // assumes allocating all of shared
	//mexPrintf("BlockSize %d = (int) floor((float)(%d-1024)/sizeof(%d)*%d*%d\n",BlockSize,deviceProp.sharedMemPerBlock,sizeof(float),sz,sz);
	//BlockSize = max(4, BlockSize);
	//BlockSize = min(BSZ, BlockSize);
	//int BlockSize = BSZ; //for no shared scenario, we don't care!
	

	//int Nfitspad= BlockSize*(int) ceil( (float) Nfitraw/BlockSize);
	//const int Nfitspad=(int) Nfitraw; //kernel auto aborts when thread > Nfitspad?

	//fun with timing
	LARGE_INTEGER freq;
	double togpu=0,fit=0,fromgpu=0,cleanup=0;
	LARGE_INTEGER start,stop;
	QueryPerformanceFrequency( &freq );
	QueryPerformanceCounter(&start);

	//get variables
	info.data=(float *) mxGetData(prhs[0]);
	info.PSFSigma=(float)mxGetScalar(prhs[1]); //matlab-dip_image convention
	info.iterations=(int) mxGetScalar(prhs[2]);
	info.fittype = (int) mxGetScalar(prhs[3]);

	if (nrhs>4)
		info.Ax = (float) mxGetScalar(prhs[4]);
	if (nrhs>5)
		info.Ay = (float) mxGetScalar(prhs[5]);
	if (nrhs>6)
		info.Bx = (float) mxGetScalar(prhs[6]);
	if (nrhs>7)
		info.By = (float) mxGetScalar(prhs[7]);
	if (nrhs>8)
		info.gamma = (float) mxGetScalar(prhs[8]);
	if (nrhs>9)
		info.d = (float) mxGetScalar(prhs[9]);
	if (nrhs>10)
		info.PSFSigma_y = (float) mxGetScalar(prhs[10]);
	else
		info.PSFSigma_y =info.PSFSigma;

#ifdef _DEBUG
	mexPrintf("subregion size: %d fittype: %d\n", info.sz, info.fittype);
	//mexPrintf("Number of fits: %d, Number padded fits: %d\n", Nfitraw, Nfitspad);
	mexPrintf("Number of fits: %d no padding\n", info.Nfitraw);
	//mexPrintf("Memory required: %dMB Memory available: %dMB\n",requiredMemory/(1024*1024),availableMemory);
#endif

	//create output for parameters and CRLBs
	switch(info.fittype){
	case 1: // (x,y,bg,I)
		plhs[0]=mxCreateNumericMatrix(info.Nfitraw, NV_P, mxSINGLE_CLASS, mxREAL);
		plhs[1]=mxCreateNumericMatrix(info.Nfitraw, NV_P, mxSINGLE_CLASS, mxREAL);
		break;
	case 2: // (x,y,bg,I,Sigma)
	case 3: // (x,y,bg,I,z)
		plhs[0]=mxCreateNumericMatrix(info.Nfitraw, NV_PS, mxSINGLE_CLASS, mxREAL);
		plhs[1]=mxCreateNumericMatrix(info.Nfitraw, NV_PS, mxSINGLE_CLASS, mxREAL);
		break;      
	case 4: // (x,y,bg,I,Sx,Sy)
		plhs[0]=mxCreateNumericMatrix(info.Nfitraw, NV_PS2, mxSINGLE_CLASS, mxREAL);
		plhs[1]=mxCreateNumericMatrix(info.Nfitraw, NV_PS2, mxSINGLE_CLASS, mxREAL);
		break;
	}
	plhs[2]=mxCreateNumericMatrix(info.Nfitraw, 1, mxSINGLE_CLASS, mxREAL);
	info.Parameters = (float*) mxGetData(plhs[0]);
	info.CRLBs = (float*) mxGetData(plhs[1]);
	info.LogLikelihood = (float*) mxGetData(plhs[2]);

	//check allocations
	if (plhs[0] == NULL || plhs[1] == NULL || plhs[2] == NULL)
		mexErrMsgIdAndTxt("cGaussMLEv2:NotEnoughMemory","Could not allocate memory for results.\n");

	QueryPerformanceCounter(&stop);
	togpu = (double) (stop.QuadPart-start.QuadPart)/freq.QuadPart;
	QueryPerformanceCounter(&start);

		SYSTEM_INFO sysinfo;
	GetSystemInfo( &sysinfo );

	const int numCPU = sysinfo.dwNumberOfProcessors;
	const int threadlimit = std::thread::hardware_concurrency();

#ifdef _DEBUG
	mexPrintf("cores %d thread limit %d\n", numCPU,threadlimit);
#endif

	std::thread *t = new std::thread[numCPU];
	for (size_t ii=0; ii<numCPU; ii++)
		t[ii] = std::thread(worker, ii, std::ref(info));
	
	for (size_t ii=0; ii<numCPU; ii++)
		t[ii].join();

	QueryPerformanceCounter(&stop);    
	fit = (double)(stop.QuadPart-start.QuadPart)/freq.QuadPart;
	QueryPerformanceCounter(&start);

	QueryPerformanceCounter(&stop);
	fromgpu = (double)(stop.QuadPart-start.QuadPart)/freq.QuadPart;
	QueryPerformanceCounter(&start);

	//cleanup
	QueryPerformanceCounter(&stop);
#ifdef _DEBUG
	//mexPrintf("Memory copies to GPU %f seconds\n", togpu);
	mexPrintf("Actual fitting time %f Actual fits per second %f\n", fit,
		Nfitraw/fit);
	//mexPrintf("Memory copies from GPU %f seconds\n", fromgpu);
	//mexPrintf("Clean up from GPU %f seconds\n", (double)(stop.QuadPart-start.QuadPart)/freq.QuadPart);
#endif
	return;
}
