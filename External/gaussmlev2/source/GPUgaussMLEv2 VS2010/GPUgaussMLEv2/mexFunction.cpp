/*!
 * \file mexFunction.cpp
 * \author Keith Lidke
 * \date January 10, 2010
 * \brief This is a pure C file that contains only the mexFunction call and the wrappers 
 * necessary to make Cuda calls.  
 * [Parameters CRLBs LL]=GPUgaussMLE(data,PSFSigma,iterations, [fit type], [Ax], [Ay], [Bx], [By], [gamma], [d], [PSFSigma_x], [PSFSigma_y])
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

extern "C" void kernel_MLEFit_wrapper(dim3 dimGrid, dim3 dimBlock, float *d_data, float PSFSigma, int sz, int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,int Nfits); 

extern "C" void kernel_MLEFit_sigma_wrapper(dim3 dimGrid, dim3 dimBlock, float *d_data, float PSFSigma, int sz, int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,int Nfits);

extern "C" void kernel_MLEFit_z_wrapper(dim3 dimGrid, dim3 dimBlock, float *d_data, float PSFSigma_x, float Ax, float Ay, float Bx, 
		float By, float gamma, float d, float PSFSigma_y, int sz, int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,int Nfits);

extern "C" void kernel_MLEFit_sigmaxy_wrapper(dim3 dimGrid, dim3 dimBlock, float *d_data, float PSFSigma, int sz, int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,int Nfits); 

//*******************************************************************************************
void CUDAERRROR(const char *instr) {
/*!
 *  \brief A simple function to dump Cuda errors and abort the run.
 *  \param instr unused string pointer
 */
    cudaError_t errornum;
    const char *str;
    if (errornum=cudaGetLastError()) {
        str=cudaGetErrorString(errornum);
        mexPrintf("%s\n", str);
        mexPrintf("You should clear this function in MATLAB for proper operation.\n", str);
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

    int blockx;
    int threadx;
    const mwSize *datasize;
    float PSFSigma, Ax, Ay, Bx, By, gamma, d, PSFSigma_y;
    int iterations, Nfits;
    float *data, *d_data;
    float *d_Parameters,*d_CRLBs,*d_LogLikelihood;
    size_t Ndim, Nfitraw, fittype, sz;

    //input checks
    if (nrhs<3)
        mexErrMsgTxt("Data input must include: data,PSFSigma,iterations,fittype !\n");
    
    if (mxGetClassID(prhs[0])!=mxSINGLE_CLASS)
        mexErrMsgTxt("Data must be comprised of single floats!\n");
    
    datasize=mxGetDimensions(prhs[0]);
    Ndim=mxGetNumberOfDimensions(prhs[0]);
    
    if (Ndim==2)Nfitraw=1;else Nfitraw=datasize[2];
    
    if (datasize[0] > IMSZBIG)
        mexErrMsgTxt("X,Y dimension of data must be smaller than 21.\n");
    if (datasize[1]!=datasize[0])
        mexErrMsgTxt("Fit Box must be square");
    sz=datasize[0];
    // mexPrintf("c sizeX: %d Ndims: %d\n",sz,Ndim);
    
    //put as many images as fit into a block
    //float szsqaure = sz*sz;
    int BlockSize = (int) floor((float)15000/4/sz/sz);
    BlockSize = max(4, BlockSize);
    BlockSize = min(BSZ, BlockSize);
    
    //mexPrintf("cdasdas %d\n",BlockSize);
    
    Nfits= BlockSize*(int) ceil( (float) Nfitraw/BlockSize);
    Nfits=(int) Nfitraw;
    //mexPrintf("c sizeX: %d Nfits: %d\n", sz, Nfits);
    //mexPrintf("Nfits: %d, Nfitraw: %d\n", Nfits, Nfitraw);

	//fun with timing
	LARGE_INTEGER freq;
	double togpu,fit,fromgpu,cleanup;
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
    
    //create device variable for data and copy to device
    cudaMalloc((void**)&d_data, sz*sz*Nfits*sizeof(float));
    cudaMemset(d_data, 0, sz*sz*Nfits*sizeof(float));
    cudaMemcpy(d_data, data, sz*sz*Nfitraw*sizeof(float), cudaMemcpyHostToDevice);
    
    //create output for parameters and CRLBs
    switch(fittype){
        case 1: // (x,y,bg,I)
            cudaMalloc((void**)&d_Parameters,   NV_P*Nfits*sizeof(float));
            cudaMalloc((void**)&d_CRLBs,        NV_P*Nfits*sizeof(float));
            plhs[0]=mxCreateNumericMatrix(Nfitraw, NV_P, mxSINGLE_CLASS, mxREAL);
            plhs[1]=mxCreateNumericMatrix(Nfitraw, NV_P, mxSINGLE_CLASS, mxREAL);
            break;
        case 2: // (x,y,bg,I,Sigma)
        case 3: // (x,y,bg,I,z)
            cudaMalloc((void**)&d_Parameters,   NV_PS*Nfits*sizeof(float));
            cudaMalloc((void**)&d_CRLBs,        NV_PS*Nfits*sizeof(float));
            plhs[0]=mxCreateNumericMatrix(Nfitraw, NV_PS, mxSINGLE_CLASS, mxREAL);
            plhs[1]=mxCreateNumericMatrix(Nfitraw, NV_PS, mxSINGLE_CLASS, mxREAL);
            break;      
        case 4: // (x,y,bg,I,Sx,Sy)
            cudaMalloc((void**)&d_Parameters,   NV_PS2*Nfits*sizeof(float));
            cudaMalloc((void**)&d_CRLBs,        NV_PS2*Nfits*sizeof(float));
            plhs[0]=mxCreateNumericMatrix(Nfitraw, NV_PS2, mxSINGLE_CLASS, mxREAL);
            plhs[1]=mxCreateNumericMatrix(Nfitraw, NV_PS2, mxSINGLE_CLASS, mxREAL);
            break;
    }
    cudaMalloc((void**)&d_LogLikelihood,        Nfits*sizeof(float));
    plhs[2]=mxCreateNumericMatrix(Nfitraw, 1, mxSINGLE_CLASS, mxREAL);
    
	QueryPerformanceCounter(&stop);
	togpu = (double) (stop.QuadPart-start.QuadPart)/freq.QuadPart;
	QueryPerformanceCounter(&start);

    //setup kernel
    blockx = (int) ceil( (float)Nfits/(float)BlockSize);
    threadx= BlockSize;
    
    dim3 dimBlock(threadx);
    dim3 dimGrid(blockx);
    
    //printf("threadx: %d,blockx: %d,Nfitraw: %d\n", threadx, blockx, Nfitraw);
    
    switch(fittype) {
        case 1: //fit x,y,bg,I
            kernel_MLEFit_wrapper(dimGrid, dimBlock, d_data, PSFSigma, (int) sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits);
            break;
            
        case 2: //fit x,y,bg,I,sigma
            kernel_MLEFit_sigma_wrapper(dimGrid, dimBlock, d_data, PSFSigma, (int) sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits);
            break;
            
        case 3: //fit x,y,bg,I,z
            kernel_MLEFit_z_wrapper(dimGrid, dimBlock, d_data, PSFSigma, Ax, Ay, Bx, By, gamma, d, PSFSigma_y, (int) sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits);
            //mexPrintf("A:%f B:%f gamma:%f d:%f \n",Ax,Bx,gamma,d);
            break;
            
        case 4: //fit x,y,bg,I,sigmax,sigmay
            kernel_MLEFit_sigmaxy_wrapper(dimGrid, dimBlock, d_data, PSFSigma, (int) sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits);
            
    }
    CUDAERRROR("kernel");
    cudaThreadSynchronize();
	QueryPerformanceCounter(&stop);    
	fit = (double)(stop.QuadPart-start.QuadPart)/freq.QuadPart;
//		Nfitraw/((double)(stop.QuadPart-start.QuadPart)/freq.QuadPart);
	QueryPerformanceCounter(&start);

	//copy to matlab output
    switch(fittype){
        case 1: // (x,y,bg,I)
            cudaMemcpy((float *)mxGetData(plhs[0]), d_Parameters, NV_P*Nfitraw*sizeof(float), cudaMemcpyDeviceToHost);
            cudaMemcpy((float *)mxGetData(plhs[1]), d_CRLBs, NV_P*Nfitraw*sizeof(float), cudaMemcpyDeviceToHost);
            break;
        case 2: // (x,y,bg,I,Sigma)
        case 3: // (x,y,bg,I,z)
            cudaMemcpy((float *)mxGetData(plhs[0]), d_Parameters, NV_PS*Nfitraw*sizeof(float), cudaMemcpyDeviceToHost);
            cudaMemcpy((float *)mxGetData(plhs[1]), d_CRLBs, NV_PS*Nfitraw*sizeof(float), cudaMemcpyDeviceToHost);
            break;
        case 4: // (x,y,bg,I,Sx,Sy)
            cudaMemcpy((float *)mxGetData(plhs[0]), d_Parameters, NV_PS2*Nfitraw*sizeof(float), cudaMemcpyDeviceToHost);
            cudaMemcpy((float *)mxGetData(plhs[1]), d_CRLBs, NV_PS2*Nfitraw*sizeof(float), cudaMemcpyDeviceToHost);
            break;
    }
    cudaMemcpy((float *)mxGetData(plhs[2]), d_LogLikelihood, Nfitraw*sizeof(float), cudaMemcpyDeviceToHost);

	QueryPerformanceCounter(&stop);
	fromgpu = (double)(stop.QuadPart-start.QuadPart)/freq.QuadPart;
	QueryPerformanceCounter(&start);
    
    //cleanup
    cudaFree(d_Parameters);
    cudaFree(d_CRLBs);
    cudaFree(d_LogLikelihood);

	cudaThreadExit(); //release context so future cudaSetDevice calls work
	QueryPerformanceCounter(&stop);
	//mexPrintf("Memory copies to GPU %f seconds\n", togpu);
    //mexPrintf("Actual fitting time %f Actual fits per second %f\n", fit,
	//	Nfitraw/fit);
	//mexPrintf("Memory copies from GPU %f seconds\n", fromgpu);
	//mexPrintf("Clean up from GPU %f seconds\n", (double)(stop.QuadPart-start.QuadPart)/freq.QuadPart);

	return;
}
