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
#ifdef _WIN32
    #include <windows.h>
#endif

#pragma comment(lib, "kernel32.lib")

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>
#include "definitions.h"
#include "CPUmleFit_LM.h"

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
	size_t Ndim, Nfitraw, fittype, sz, cameratype;

	//sCMOS
	float *varim=0, *d_varim=0, *gainim=0, *d_gainim=0;

	//Spline
	int spline_xsize, spline_ysize, spline_zsize;
	float *coeff=0,*d_coeff=0;
	const mwSize *datasize_spline;

	int silent = 0;

	//input checks
	if (nrhs<3)
		mexErrMsgIdAndTxt("cGaussMLEv2:WrongNoArgs","Inputs must include: cGaussMLEv2(data,PSFSigma,iterations,fittype) !\n");

	if (mxGetClassID(prhs[0])!=mxSINGLE_CLASS)
		mexErrMsgIdAndTxt("cGaussMLEv2:WrongDataType","Data must be comprised of single floats!\n");

	silent = (int) mxGetScalar(prhs[5]);
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
	mwSize freq;
	double togpu=0,fit=0,fromgpu=0,cleanup=0;
	mwSize start,stop;
    #ifdef _WIN32
        QueryPerformanceFrequency( &freq );
        QueryPerformanceCounter(&start);
    #endif
    


	//get variables
	data=(float *) mxGetData(prhs[0]);
	//PSFSigma=(float)mxGetScalar(prhs[1]); //matlab-dip_image convention
	iterations=(int) mxGetScalar(prhs[2]);
	fittype = (int) mxGetScalar(prhs[3]);
	cameratype = (int) mxGetScalar(prhs[4]);

	if (fittype != 5){
		PSFSigma=(float)mxGetScalar(prhs[1]);
		if (cameratype ==0){//EMCCD
			if (nrhs>6)
				Ax = (float) mxGetScalar(prhs[6]);
			if (nrhs>7)
				Ay = (float) mxGetScalar(prhs[7]);
			if (nrhs>8)
				Bx = (float) mxGetScalar(prhs[8]);
			if (nrhs>9)
				By = (float) mxGetScalar(prhs[9]);
			if (nrhs>10)
				gamma = (float) mxGetScalar(prhs[10]);
			if (nrhs>11)
				d = (float) mxGetScalar(prhs[11]);
			if (nrhs>12)
				PSFSigma_y = (float) mxGetScalar(prhs[12]);
			else
				PSFSigma_y =PSFSigma;
		}
		else if (cameratype ==1){//sCMOS
			varim=(float *) mxGetData(prhs[6]);
			gainim=(float *) mxGetData(prhs[7]);
			if (nrhs>8)
				Ax = (float) mxGetScalar(prhs[8]);
			if (nrhs>9)
				Ay = (float) mxGetScalar(prhs[9]);
			if (nrhs>10)
				Bx = (float) mxGetScalar(prhs[10]);
			if (nrhs>11)
				By = (float) mxGetScalar(prhs[11]);
			if (nrhs>12)
				gamma = (float) mxGetScalar(prhs[12]);
			if (nrhs>13)
				d = (float) mxGetScalar(prhs[13]);
			if (nrhs>14)
				PSFSigma_y = (float) mxGetScalar(prhs[14]);
		}
	} 

	if (fittype == 5){
		coeff =(float *) mxGetData(prhs[1]);
		datasize_spline=mxGetDimensions(prhs[1]);
		spline_xsize = datasize_spline[0];
		spline_ysize = datasize_spline[1];
		spline_zsize = datasize_spline[2];
		if (cameratype ==1){//sCMOS
			varim=(float *) mxGetData(prhs[6]);
			gainim=(float *) mxGetData(prhs[7]);
		}
		//cudaMalloc((void**)&d_coeff, spline_xsize*spline_ysize*spline_zsize*64*sizeof(float));
		//cudaMemset(d_coeff, 0, spline_xsize*spline_ysize*spline_zsize*64*sizeof(float));
		//cudaMemcpy(d_coeff, coeff, spline_xsize*spline_ysize*spline_zsize*64*sizeof(float), cudaMemcpyHostToDevice);
	}
if  (silent==0)
{
	mexPrintf("subregion size: %d fittype: %d\n", sz, fittype);
	//mexPrintf("Number of fits: %d, Number padded fits: %d\n", Nfitraw, Nfitspad);
	mexPrintf("Number of fits: %d no padding\n", Nfitraw);
	//mexPrintf("Memory required: %dMB Memory available: %dMB\n",requiredMemory/(1024*1024),availableMemory);
}

	//create output for parameters and CRLBs
	switch(fittype){
	case 1: // (x,y,bg,I)
		plhs[0]=mxCreateNumericMatrix(Nfitraw, NV_P+1, mxSINGLE_CLASS, mxREAL);
		plhs[1]=mxCreateNumericMatrix(Nfitraw, NV_P, mxSINGLE_CLASS, mxREAL);
		break;
	case 2: // (x,y,bg,I,Sigma)
	case 3: // (x,y,bg,I,z)
		plhs[0]=mxCreateNumericMatrix(Nfitraw, NV_PS+1, mxSINGLE_CLASS, mxREAL);
		plhs[1]=mxCreateNumericMatrix(Nfitraw, NV_PS, mxSINGLE_CLASS, mxREAL);
		break;      
	case 4: // (x,y,bg,I,Sx,Sy)
		plhs[0]=mxCreateNumericMatrix(Nfitraw, NV_PS2+1, mxSINGLE_CLASS, mxREAL);
		plhs[1]=mxCreateNumericMatrix(Nfitraw, NV_PS2, mxSINGLE_CLASS, mxREAL);
		break;
	case 5:
		plhs[0]=mxCreateNumericMatrix(Nfitraw, NV_PSP+1, mxSINGLE_CLASS, mxREAL);
		plhs[1]=mxCreateNumericMatrix(Nfitraw, NV_PSP, mxSINGLE_CLASS, mxREAL);
		break;
	}
	plhs[2]=mxCreateNumericMatrix(Nfitraw, 1, mxSINGLE_CLASS, mxREAL);
	Parameters = (float*) mxGetData(plhs[0]);
	CRLBs = (float*) mxGetData(plhs[1]);
	LogLikelihood = (float*) mxGetData(plhs[2]);

	//check allocations
	if (plhs[0] == NULL || plhs[1] == NULL || plhs[2] == NULL)
		mexErrMsgIdAndTxt("cGaussMLEv2:NotEnoughMemory","Could not allocate memory for results.\n");
    #ifdef _WIN32
        QueryPerformanceCounter(&stop);
        togpu = (double) (stop.QuadPart-start.QuadPart)/freq.QuadPart;
        QueryPerformanceCounter(&start);
    #endif


	//setup kernel
//	blockx = (int) ceil( (float)Nfitraw/(float)BlockSize);
//	threadx= BlockSize;
//
//#ifdef _DEBUG
//	mexPrintf("threadx: %d,blockx: %d\n", threadx, blockx);
//#endif

	switch(fittype) {
	case 1: //fit x,y,bg,I
		for (int ii = 0; ii<Nfitraw; ii++)
			kernel_MLEFit_LM_EMCCD(ii, data, PSFSigma, (int) sz, iterations, Parameters, CRLBs, LogLikelihood, (const int) Nfitraw);
		//kernel_MLEFit_wrapper(dimGrid, dimBlock, d_data, PSFSigma, (int) sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfitraw);
		break;

	case 2: //fit x,y,bg,I,sigma
		for (int ii = 0; ii<Nfitraw; ii++)
			kernel_MLEFit_LM_Sigma_EMCCD(ii, data, PSFSigma, (int) sz, iterations, Parameters, CRLBs, LogLikelihood, (const int) Nfitraw);
		break;

	case 3: //fit x,y,bg,I,z
		for (int ii = 0; ii<Nfitraw; ii++)
			kernel_MLEFit_LM_z_EMCCD(ii, data, PSFSigma, Ax, Ay, Bx, By, gamma, d, PSFSigma_y, (int) sz, iterations, Parameters, CRLBs, LogLikelihood, (const int) Nfitraw);
		//mexPrintf("A:%f B:%f gamma:%f d:%f \n",Ax,Bx,gamma,d);
		break;

	case 4: //fit x,y,bg,I,sigmax,sigmay
		for (int ii = 0; ii<Nfitraw; ii++)
			kernel_MLEFit_LM_sigmaxy_EMCCD(ii, data, PSFSigma, (int) sz, iterations, Parameters, CRLBs, LogLikelihood, (const int) Nfitraw);
		break;
	case 5: //fit x,y,bg,I,sigmax,sigmay
		for (int ii = 0; ii<Nfitraw; ii++)
			kernel_splineMLEFit_z_EMCCD(ii,data,coeff,spline_xsize,spline_ysize,spline_zsize,(int) sz,iterations,Parameters,CRLBs,LogLikelihood,(const int) Nfitraw);

			//kernel_MLEFit_LM_sigmaxy_EMCCD(ii, data, PSFSigma, (int) sz, iterations, Parameters, CRLBs, LogLikelihood, (const int) Nfitraw);

	}
    #ifdef _WIN32
	QueryPerformanceCounter(&stop);    
	fit = (double)(stop.QuadPart-start.QuadPart)/freq.QuadPart;
	QueryPerformanceCounter(&start);

	QueryPerformanceCounter(&stop);
	fromgpu = (double)(stop.QuadPart-start.QuadPart)/freq.QuadPart;
	QueryPerformanceCounter(&start);

	//cleanup
	QueryPerformanceCounter(&stop);
    #endif  

	if (silent==0)
	{//mexPrintf("Memory copies to GPU %f seconds\n", togpu);
		mexPrintf("Actual fitting time %f Actual fits per second %f\n", fit,
			Nfitraw/fit);
		//mexPrintf("Memory copies from GPU %f seconds\n", fromgpu);
		//mexPrintf("Clean up from GPU %f seconds\n", (double)(stop.QuadPart-start.QuadPart)/freq.QuadPart);
	}
	return;
}
