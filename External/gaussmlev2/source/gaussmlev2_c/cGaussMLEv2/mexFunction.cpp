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
	float *data=0;
	float *Parameters=0,*CRLBs=0,*LogLikelihood=0;
	size_t Ndim, Nfitraw, fittype, sz;

	//input checks
	if (nrhs<3)
		mexErrMsgIdAndTxt("cGaussMLEv2:WrongNoArgs","Inputs must include: cGaussMLEv2(data,PSFSigma,iterations,fittype) !\n");

	if (mxGetClassID(prhs[0])!=mxSINGLE_CLASS)
		mexErrMsgIdAndTxt("cGaussMLEv2:WrongDataType","Data must be comprised of single floats!\n");

	datasize=mxGetDimensions(prhs[0]);
	Ndim=mxGetNumberOfDimensions(prhs[0]);

	if (Ndim==2)Nfitraw=1;else Nfitraw=datasize[2];

	if (datasize[0] > IMSZBIG)
		mexErrMsgIdAndTxt("cGaussMLEv2:SubregionSize","X,Y dimension of data must be smaller than 21.\n");
	if (datasize[1]!=datasize[0])
		mexErrMsgIdAndTxt("cGaussMLEv2:SubregionShape","Fit Box must be square");
	sz=datasize[0];
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

#ifdef _DEBUG
	mexPrintf("subregion size: %d fittype: %d\n", sz, fittype);
	//mexPrintf("Number of fits: %d, Number padded fits: %d\n", Nfitraw, Nfitspad);
	mexPrintf("Number of fits: %d no padding\n", Nfitraw);
	//mexPrintf("Memory required: %dMB Memory available: %dMB\n",requiredMemory/(1024*1024),availableMemory);
#endif

	//create output for parameters and CRLBs
	switch(fittype){
	case 1: // (x,y,bg,I)
		plhs[0]=mxCreateNumericMatrix(Nfitraw, NV_P, mxSINGLE_CLASS, mxREAL);
		plhs[1]=mxCreateNumericMatrix(Nfitraw, NV_P, mxSINGLE_CLASS, mxREAL);
		break;
	case 2: // (x,y,bg,I,Sigma)
	case 3: // (x,y,bg,I,z)
		plhs[0]=mxCreateNumericMatrix(Nfitraw, NV_PS, mxSINGLE_CLASS, mxREAL);
		plhs[1]=mxCreateNumericMatrix(Nfitraw, NV_PS, mxSINGLE_CLASS, mxREAL);
		break;      
	case 4: // (x,y,bg,I,Sx,Sy)
		plhs[0]=mxCreateNumericMatrix(Nfitraw, NV_PS2, mxSINGLE_CLASS, mxREAL);
		plhs[1]=mxCreateNumericMatrix(Nfitraw, NV_PS2, mxSINGLE_CLASS, mxREAL);
		break;
	}
	plhs[2]=mxCreateNumericMatrix(Nfitraw, 1, mxSINGLE_CLASS, mxREAL);
	Parameters = (float*) mxGetData(plhs[0]);
	CRLBs = (float*) mxGetData(plhs[1]);
	LogLikelihood = (float*) mxGetData(plhs[2]);

	//check allocations
	if (plhs[0] == NULL || plhs[1] == NULL || plhs[2] == NULL)
		mexErrMsgIdAndTxt("cGaussMLEv2:NotEnoughMemory","Could not allocate memory for results.\n");

	QueryPerformanceCounter(&stop);
	togpu = (double) (stop.QuadPart-start.QuadPart)/freq.QuadPart;
	QueryPerformanceCounter(&start);

	//setup kernel
	blockx = (int) ceil( (float)Nfitraw/(float)BlockSize);
	threadx= BlockSize;

#ifdef _DEBUG
	mexPrintf("threadx: %d,blockx: %d\n", threadx, blockx);
#endif

	switch(fittype) {
	case 1: //fit x,y,bg,I
		for (int ii = 0; ii<Nfitraw; ii++)
			kernel_MLEFit(ii, data, PSFSigma, (int) sz, iterations, Parameters, CRLBs, LogLikelihood, (const int) Nfitraw);
		//kernel_MLEFit_wrapper(dimGrid, dimBlock, d_data, PSFSigma, (int) sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfitraw);
		break;

	case 2: //fit x,y,bg,I,sigma
		for (int ii = 0; ii<Nfitraw; ii++)
			kernel_MLEFit_sigma(ii, data, PSFSigma, (int) sz, iterations, Parameters, CRLBs, LogLikelihood, (const int) Nfitraw);
		break;

	case 3: //fit x,y,bg,I,z
		for (int ii = 0; ii<Nfitraw; ii++)
			kernel_MLEFit_z(ii, data, PSFSigma, Ax, Ay, Bx, By, gamma, d, PSFSigma_y, (int) sz, iterations, Parameters, CRLBs, LogLikelihood, (const int) Nfitraw);
		//mexPrintf("A:%f B:%f gamma:%f d:%f \n",Ax,Bx,gamma,d);
		break;

	case 4: //fit x,y,bg,I,sigmax,sigmay
		for (int ii = 0; ii<Nfitraw; ii++)
			kernel_MLEFit_sigmaxy(ii, data, PSFSigma, (int) sz, iterations, Parameters, CRLBs, LogLikelihood, (const int) Nfitraw);

	}
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
