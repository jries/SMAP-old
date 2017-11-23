//Copyright (c) 2017 Ries Lab, European Molecular Biology Laboratory, Heidelberg.
//author: Yiming Li
//email: yiming.li@embl.de
//date: 2017.07.24

//Using Keith Lidke's template
//"Fast, single-molecule localization that achieves theoretically minimum uncertainty." 
//C. Simth, N. Joseph, B. Rieger & K. Lidke. Nat. Methods, 7, 373, 2010

/*! 
* \file mexFunction.cpp
* \brief This is a pure C file that contains only the mexFunction call and the wrappers 
* necessary to make Cuda calls.  
* [Parameters CRLBs LL]=CPUmleFit(varargin)
*  \varargin:
*  \imstack, startpsf/coeff, iterations, fitmode, isemccd, hidereport
*  \1. imagestack (single)
*  \2. fitmode
*  \1 fix PSF
*  \2 free PSF
*  \3 Gauss fit z
*  \4 fit PSFx, PSFy elliptical
*  \5 cspline
*  \optional:
*  \3. iterations (default=50)
*  \4. paramters for fitters:
*  \1. fix PSF: PSFxy sigma
*  \2. free PSF: start PSFxy
*  \3. Gauss fit z: parameters: PSFx, Ax, Ay, Bx, By, gamma, d, PSFy
*  \(single)
*  \4. fit PSFx, PSFy elliptical: start PSFx, PSFy
*  \5. cspline: cspline coefficients (single)
*  \6. cspline for 2D PSF: as 5, but two fits to break asymmetry. cspline coefficients (single)
*  \5. varmap: Variance map for sCMOS. If6 size(varmap) ~= size(imagestack):
*  \no sCMOS correction is used. Default= emCCD
*  \6. silent (suppress output if 1)

*  \Output:
*  \P
*  \1. X, Y, Photons, Background, Iterations
*  \2. X, Y, Photons, Background, PSFxy, Iterations
*  \3. X, Y, Photons, Background, Z, Iterations
*  \4. X, Y, Photons, Background, PSFx, PSFy, Iterations
*  \5. X, Y, Photons, Background, Z, Iterations
*  \6. X, Y, Photons, Background, Z, Iterations
*  \CRLB: cramer-rao lower bounds, as in P
*  \LogL: log-likelihood.

//Terms of Use 
//
//This file is part of GPUmleFit_LM. 
//
//GPUmleFit_LM Fitter is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. 
//
//GPUmleFit_LM Fitter is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. 
//
//You should have received a copy of the GNU General Public License along with GPUmleFit_LM Fitter. If not, see <http://www.gnu.org/licenses/>. 
//
//Additional permission under GNU GPL version 3 section 7 

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
    float * test;
	int blockx=0;
	int threadx=0;
	const mwSize *datasize=0;
	//const int *datasize=0;
	int iterations=50,i;
	float * startParameters = 0;
	float *data=0, *d_data=0, *dTAll, *d_dTAll;
	float *Parameters=0,*CRLBs=0,*LogLikelihood=0, *residual=0;
	//size_t Ndim, fittype, cameratype=0,noParameters,Nfitraw,noChannels,sz;;//default cameratype 0 EMCCD
	int Ndim, fittype, cameratype=0,noParameters,Nfitraw,noChannels,sz;;//default cameratype 0 EMCCD
	//int noParameters,Nfitraw,noChannels,sz;
	int sumShared=0;
	int *shared, *d_shared;
	//sCMOS
	float *varim=0, *d_varim=0;

	//Spline
	int spline_xsize, spline_ysize, spline_zsize;
	float *coeff=0,*d_coeff=0;
	const mwSize *datasize_spline;
	//const int *datasize_spline;
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
	//if (nrhs<2)
	//	mexErrMsgIdAndTxt("CPUmleFit_LM:WrongNoArgs","Inputs must include: CPUmleFit_LM(data,fittype,iterations,startparameter,var,silent) !\n");
	//else
	//{
	if (mxGetClassID(prhs[0])!=mxSINGLE_CLASS)
		mexErrMsgIdAndTxt("CPUmleFit_LM:WrongDataType","Data must be comprised of single floats!\n");

	datasize=mxGetDimensions(prhs[0]);
	Ndim=mxGetNumberOfDimensions(prhs[0]);

	Nfitraw=datasize[2];

	if (datasize[0] > IMSZBIG)
		mexErrMsgIdAndTxt("CPUmleFit_LM:SubregionSize","X,Y dimension of data must be smaller than 51.\n");
	if (datasize[1]!=datasize[0])
		mexErrMsgIdAndTxt("CPUmleFit_LM:SubregionShape","Fit Box must be square");
	sz=datasize[0];
	noChannels = datasize[3];

	
	data=(float *) mxGetData(prhs[0]);

	shared = (int *) mxGetData(prhs[1]);

	iterations=(int) mxGetScalar(prhs[2]);



	if (mxGetClassID(prhs[3])!=mxSINGLE_CLASS)
		mexErrMsgIdAndTxt("GPUmleFit_LM:WrongDataType","coeff must be comprised of single floats!\n");

	coeff =(float *) mxGetData(prhs[3]);
	datasize_spline=mxGetDimensions(prhs[3]);
	spline_xsize = datasize_spline[0];
	spline_ysize = datasize_spline[1];
	spline_zsize = datasize_spline[2];

	dTAll = (float *) mxGetData(prhs[4]);


	sumShared=0;
	for (i=0;i<5;i++){
		sumShared+=shared[i];
	}
	noParameters = (5*noChannels-sumShared*(noChannels-1));

	/*mexPrintf("noChannels %f noParameters %f sumShared %f Nfitraw%f test%f \n", noChannels,noParameters,sumShared, Nfitraw, 1);*/


	plhs[0]=mxCreateNumericMatrix(Nfitraw, noParameters+1, mxSINGLE_CLASS, mxREAL);
	plhs[1]=mxCreateNumericMatrix(Nfitraw, noParameters, mxSINGLE_CLASS, mxREAL);

	plhs[2]=mxCreateNumericMatrix(Nfitraw, 1, mxSINGLE_CLASS, mxREAL);

	plhs[3]=mxCreateNumericArray(Ndim, datasize, mxSINGLE_CLASS, mxREAL);

	//****************************test*****************************
	//plhs[3]=mxCreateNumericMatrix(100, 1, mxSINGLE_CLASS, mxREAL);
	//test = (float*) mxGetData(plhs[3]);
	//test[0]=1;
	//test[1]=noChannels;
	//test[2]=noParameters;
	//test[3]=sumShared;
	//test[4]=Nfitraw;
	//test[5]=shared[0];
	//test[6]=shared[1];
	//test[7]=shared[2];
	//test[8]=shared[3];
	//test[9]=shared[4];
	//****************************test*****************************

	
	Parameters = (float*) mxGetData(plhs[0]);
	CRLBs = (float*) mxGetData(plhs[1]);
	LogLikelihood = (float*) mxGetData(plhs[2]);
	residual = (float*) mxGetData(plhs[3]);

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
	/*mexPrintf("before iterations %f\n",(double)iterations);*/
	for (int ii = 0; ii<Nfitraw; ii++)
		kernel_splineMLEFit_z_EMCCD_multi(ii,data,coeff,dTAll,spline_xsize,spline_ysize,spline_zsize,(int) sz,iterations,noParameters,noChannels,Parameters,CRLBs,LogLikelihood,residual,spline_zsize/2.0f,(const int) Nfitraw,shared);
			
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
