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
	float PSFSigma=1, Ax=0, Ay=0, Bx=0, By=0, gamma=0.5, d=0.5, PSFSigma_y=1;
	int iterations=50;
	float * startParameters = 0;
	float *data=0, *d_data=0;
	float *Parameters=0,*CRLBs=0,*LogLikelihood=0;
	size_t Ndim, Nfitraw, fittype, sz, cameratype=0;//default cameratype 0 EMCCD
	mwSize noParameters;
	//sCMOS
	float *varim=0;

	//Spline
	int spline_xsize, spline_ysize, spline_zsize;
	float *coeff=0,*d_coeff=0;
	const mwSize *datasize_spline;

	int silent = 0;

	//fun with timing
	//mwSize freq;
	
	//mwSize start,stop;
    #ifdef _WIN32
        LARGE_INTEGER freq;
        LARGE_INTEGER start,stop;
        double togpu=0,fit=0,fromgpu=0,cleanup=0;
        QueryPerformanceFrequency( &freq );
        QueryPerformanceCounter(&start);
    #endif


	//input checks
	if (nrhs<2)
		mexErrMsgIdAndTxt("CPUmleFit_LM:WrongNoArgs","Inputs must include: CPUmleFit_LM(data,fittype,iterations,startparameter,var,silent) !\n");
	else
	{
		if (mxGetClassID(prhs[0])!=mxSINGLE_CLASS)
		mexErrMsgIdAndTxt("CPUmleFit_LM:WrongDataType","Data must be comprised of single floats!\n");

		datasize=mxGetDimensions(prhs[0]);
		Ndim=mxGetNumberOfDimensions(prhs[0]);

		if (Ndim==2)Nfitraw=1;else Nfitraw=datasize[2];

		if (datasize[0] > IMSZBIG)
			mexErrMsgIdAndTxt("CPUmleFit_LM:SubregionSize","X,Y dimension of data must be smaller than 21.\n");
		if (datasize[1]!=datasize[0])
			mexErrMsgIdAndTxt("CPUmleFit_LM:SubregionShape","Fit Box must be square");
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
				if (mxGetClassID(prhs[3])!=mxSINGLE_CLASS)
					mexErrMsgIdAndTxt("CPUmleFit_LM:WrongNoArgs","Parameters should be single. \n");
				else
				{
					startParameters = (float*)mxGetData(prhs[3]);
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
				}
				break;
			case 5://spline fit
			case 6:
				iterations=(int) mxGetScalar(prhs[2]);
				if (mxGetClassID(prhs[3])!=mxSINGLE_CLASS)
					mexErrMsgIdAndTxt("GPUmleFit_LM:WrongDataType","coeff must be comprised of single floats!\n");

				coeff =(float *) mxGetData(prhs[3]);
				datasize_spline=mxGetDimensions(prhs[3]);
				spline_xsize = datasize_spline[0];
				spline_ysize = datasize_spline[1];
				spline_zsize = datasize_spline[2];
				break;
			//case 6:
			//	//iterations=(int) mxGetScalar(prhs[2]);
			//	if (mxGetClassID(prhs[3])!=mxSINGLE_CLASS)
			//		mexErrMsgIdAndTxt("GPUmleFit_LM:WrongDataType","coeff must be comprised of single floats!\n");

			//	coeff =(float *) mxGetData(prhs[3]);
			//	datasize_spline=mxGetDimensions(prhs[3]);
			//	spline_xsize = datasize_spline[0];
			//	spline_ysize = datasize_spline[1];
			//	spline_zsize = datasize_spline[2];
			//	break;
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
	case 6:
		plhs[0]=mxCreateNumericMatrix(Nfitraw, (NV_PS+2), mxSINGLE_CLASS, mxREAL);
		plhs[1]=mxCreateNumericMatrix(Nfitraw, NV_PS, mxSINGLE_CLASS, mxREAL);
		plhs[3]=mxCreateNumericMatrix(Nfitraw, (NV_PS+2), mxSINGLE_CLASS, mxREAL);
		plhs[4]=mxCreateNumericMatrix(Nfitraw, NV_PS, mxSINGLE_CLASS, mxREAL);
		plhs[5]=mxCreateNumericMatrix(Nfitraw, 1, mxSINGLE_CLASS, mxREAL);
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
			if (cameratype==0){
				for (int ii = 0; ii<Nfitraw; ii++)
					kernel_MLEFit_LM_EMCCD(ii, data, PSFSigma, (int) sz, iterations, Parameters, CRLBs, LogLikelihood, (const int) Nfitraw);
			}
			else if (cameratype ==1) {
				for (int ii = 0; ii<Nfitraw; ii++)
					kernel_MLEFit_LM_sCMOS(ii, data, PSFSigma, (int) sz, iterations, Parameters, CRLBs, LogLikelihood, (const int) Nfitraw,varim);
			}
			//kernel_MLEFit_wrapper(dimGrid, dimBlock, d_data, PSFSigma, (int) sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfitraw);
			break;

		case 2: //fit x,y,bg,I,sigma
			if (cameratype==0){
				for (int ii = 0; ii<Nfitraw; ii++)
					kernel_MLEFit_LM_Sigma_EMCCD(ii, data, PSFSigma, (int) sz, iterations, Parameters, CRLBs, LogLikelihood, (const int) Nfitraw);
			}
			else if (cameratype ==1) {
				for (int ii = 0; ii<Nfitraw; ii++)
					kernel_MLEFit_LM_Sigma_sCMOS(ii, data, PSFSigma, (int) sz, iterations, Parameters, CRLBs, LogLikelihood, (const int) Nfitraw,varim);
			}
			break;

		case 3: //fit x,y,bg,I,z
			if (cameratype==0){
				for (int ii = 0; ii<Nfitraw; ii++)
					kernel_MLEFit_LM_z_EMCCD(ii, data, PSFSigma, Ax, Ay, Bx, By, gamma, d, PSFSigma_y, (int) sz, iterations, Parameters, CRLBs, LogLikelihood, (const int) Nfitraw);
			}
			else if (cameratype ==1){
				for (int ii = 0; ii<Nfitraw; ii++)
					kernel_MLEFit_LM_z_sCMOS(ii, data, PSFSigma, Ax, Ay, Bx, By, gamma, d, PSFSigma_y, (int) sz, iterations, Parameters, CRLBs, LogLikelihood, (const int) Nfitraw,varim);
			}
			//mexPrintf("A:%f B:%f gamma:%f d:%f \n",Ax,Bx,gamma,d);
			break;

		case 4: //fit x,y,bg,I,sigmax,sigmay
			if (cameratype==0){
			for (int ii = 0; ii<Nfitraw; ii++)
				kernel_MLEFit_LM_sigmaxy_EMCCD(ii, data, PSFSigma, (int) sz, iterations, Parameters, CRLBs, LogLikelihood, (const int) Nfitraw);
			}
			else if (cameratype ==1){
				for (int ii = 0; ii<Nfitraw; ii++)
					kernel_MLEFit_LM_sigmaxy_sCMOS(ii, data, PSFSigma, (int) sz, iterations, Parameters, CRLBs, LogLikelihood, (const int) Nfitraw,varim);
			}
			break;
		case 5: //fit x,y,bg,I,sigmax,sigmay
			if (cameratype==0){
				for (int ii = 0; ii<Nfitraw; ii++)
					kernel_splineMLEFit_z_EMCCD(ii,data,coeff,spline_xsize,spline_ysize,spline_zsize,(int) sz,iterations,Parameters,CRLBs,LogLikelihood,(const int) Nfitraw);
			}
			else if (cameratype ==1){
				for (int ii = 0; ii<Nfitraw; ii++)
					kernel_splineMLEFit_z_sCMOS(ii,data,coeff,spline_xsize,spline_ysize,spline_zsize,(int) sz,iterations,Parameters,CRLBs,LogLikelihood,(const int) Nfitraw,varim);
			}
			break;
		case 6:
			if (cameratype==0){
				for (int ii = 0; ii<Nfitraw; ii++)
				kernel_splineMLEFit2D_z_EMCCD(ii, data, coeff, spline_xsize, spline_ysize,  spline_zsize,  (int)sz,  iterations, 
					Parameters,  CRLBs,  LogLikelihood, spline_zsize/4.0f*3.0f,(const int)Nfitraw);

				for (int ii = 0; ii<Nfitraw; ii++)
				kernel_splineMLEFit2D_z_EMCCD(ii, data, coeff, spline_xsize, spline_ysize,  spline_zsize,  (int)sz,  iterations, 
					(float*)mxGetData(plhs[3]),  (float*)mxGetData(plhs[4]),  (float*)mxGetData(plhs[5]), spline_zsize/4.0f,(const int)Nfitraw);
			}
			else if (cameratype ==1){
				for (int ii = 0; ii<Nfitraw; ii++)
				kernel_splineMLEFit2D_z_sCMOS(ii, data, coeff, spline_xsize, spline_ysize,  spline_zsize,  (int)sz,  iterations, 
					Parameters,  CRLBs,  LogLikelihood, spline_zsize/4.0f*3.0f,(const int)Nfitraw,varim);

				for (int ii = 0; ii<Nfitraw; ii++)
				kernel_splineMLEFit2D_z_sCMOS(ii, data, coeff, spline_xsize, spline_ysize,  spline_zsize,  (int)sz,  iterations, 
					(float*)mxGetData(plhs[3]),  (float*)mxGetData(plhs[4]),  (float*)mxGetData(plhs[5]), spline_zsize/4.0f,(const int)Nfitraw,varim);
				
			}
			break;

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
     

	if (silent==0)
	{//mexPrintf("Memory copies to GPU %f seconds\n", togpu);
		mexPrintf("Actual fitting time %f Actual fits per second %f\n", fit,
			Nfitraw/fit);
		//mexPrintf("Memory copies from GPU %f seconds\n", fromgpu);
		//mexPrintf("Clean up from GPU %f seconds\n", (double)(stop.QuadPart-start.QuadPart)/freq.QuadPart);
	}
    #endif 
	return;
}
