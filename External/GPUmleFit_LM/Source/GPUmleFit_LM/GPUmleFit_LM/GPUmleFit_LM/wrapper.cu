/*!
 * \file wrapper.cu
 * \author Keith Lidke
 * \date January 10, 2010
 * \brief Wrap the Cuda kernel calls as standard external C functions.  This allows the kernels to be
 * called without doing anything special in the C code and simplifies building the code.
 */
#include "definitions.h"
#include "GPUmleFit_LM_EMCCD.h"
#include "GPUgaussLib.h"
//#include "GPUgaussLib.cuh"
//#include "GPUgaussMLE_LM.h"
//#include "GPUgaussMLE_sCMOS.h"
//#include "GPUsplineMLE.h"

//EMCCD wrapper

//*******************************************************************************************
extern void kernel_MLEFit_EMCCD_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits) 
{
/*!
 *  \brief Basic maximum likelihood estimator fit based kernel
 *  \param dimGrid number of blocks 
 *  \param dimBlock number of threads per block
 *  \param d_data an array of subregions to be processed copied into video memory
 *  \param PSFSigma the sigma value to use for the point spread function
 *  \param sz nxn size of the subregion
 *  \param iterations maximum allowed iterations before aborting fitting
 *  \param d_Parameters pointer to result array of fitted parameters, x, y coords, etc.
 *  \param d_CRLBs pointer to result array of Cramer-Rao lower bound estimates 
 *  \param d_LogLikelihood pointer to result array of loglikelihood estimates of fitting
 *  \param Nfits number of subregions to fit
 */

	kernel_MLEFit_LM_EMCCD<<<dimGrid, dimBlock>>>(d_data, PSFSigma, sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits);
	
}

//*******************************************************************************************
extern void kernel_MLEFit_sigma_EMCCD_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits) 
{
/*!
 *  \brief Basic maximum likelihood estimator fit based kernel
 *  \param dimGrid number of blocks 
 *  \param dimBlock number of threads per block
 *  \param d_data an array of subregions to be processed copied into video memory
 *  \param PSFSigma the sigma value to use for the point spread function
 *  \param sz nxn size of the subregion
 *  \param iterations maximum allowed iterations before aborting fitting
 *  \param d_Parameters pointer to result array of fitted parameters, x, y coords, etc.
 *  \param d_CRLBs pointer to result array of Cramer-Rao lower bound estimates 
 *  \param d_LogLikelihood pointer to result array of loglikelihood estimates of fitting
 *  \param Nfits number of subregions to fit
 */
	//kernel_MLEFit_sigma_EMCCD_LM<<<dimGrid, dimBlock>>>(d_data, PSFSigma, sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits);
	//void (*jacf)(int, int,float*,float*,float*) = kernel_DerivativeGauss2D_sigma;
	kernel_MLEFit_LM_Sigma_EMCCD<<<dimGrid, dimBlock>>>(d_data, PSFSigma, sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits);
	//InitPara Para;


	//Para.PSFSigma=NULL;
 //    Para.sz=NULL;
	//Para.iterations=NULL;
	// Para. Nfits=NULL;
	// Para.NV=NULL;
 //    Para.clamp=NULL;


	// Para.PSFSigma_x=NULL;
	// Para.PSFSigma_y=NULL;
	// Para.Ax=NULL;
	// Para.Ay=NULL;
	// Para. Bx=NULL;
	// Para. By=NULL;
	// Para. gamma=NULL;
	// Para. d=NULL;

	// Para. coeff=NULL;
	// Para. spline_xsize=NULL;
	// Para. spline_ysize=NULL;
	// Para. spline_zsize=NULL;




	//Para.sz = sz;
	//Para.PSFSigma=PSFSigma;
	//Para.iterations = iterations;
	//Para.Nfits = Nfits;
	//Para.NV = 5;
	//float clamp[5]={1.0,1.0,100,20,0.5};
	//Para.clamp= clamp;
	//kernel_MLEFit_LM_EMCCD<<<dimGrid, dimBlock>>>(d_data, &Para, d_Parameters, d_CRLBs, d_LogLikelihood,2);

}

//*******************************************************************************************
extern void kernel_MLEFit_z_EMCCD_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma_x, const float Ax, const float Ay, const float Bx, 
		const float By, const float gamma, const float d, const float PSFSigma_y, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits) 
{
/*!
 *  \brief Basic maximum likelihood estimator fit based kernel
 *  \param dimGrid number of blocks 
 *  \param dimBlock number of threads per block
 *  \param d_data an array of subregions to be processed copied into video memory
 *  \param PSFSigma_x the sigma value to use for the point spread function on the x axis
 *  \param Ax ???
 *  \param Ay ???
 *  \param Bx ???
 *  \param By ???
 *  \param gamma ???
 *  \param d ???
 *  \param PSFSigma_y the sigma value to use for the point spread function on the y axis
 *  \param sz nxn size of the subregion
 *  \param iterations maximum allowed iterations before aborting fitting
 *  \param d_Parameters pointer to result array of fitted parameters, x, y coords, etc.
 *  \param d_CRLBs pointer to result array of Cramer-Rao lower bound estimates 
 *  \param d_LogLikelihood pointer to result array of loglikelihood estimates of fitting
 *  \param Nfits number of subregions to fit
 */
	kernel_MLEFit_LM_z_EMCCD<<<dimGrid, dimBlock>>>(d_data, PSFSigma_x, Ax, Ay, Bx, By, gamma, d, PSFSigma_y, sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits);
}

//*******************************************************************************************
extern void kernel_MLEFit_sigmaxy_EMCCD_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits) 
{
	/*!
 *  \brief Basic maximum likelihood estimator fit based kernel
 *  \param dimGrid number of blocks 
 *  \param dimBlock number of threads per block
 *  \param d_data an array of subregions to be processed copied into video memory
 *  \param PSFSigma the sigma value to use for the point spread function
 *  \param sz nxn size of the subregion
 *  \param iterations maximum allowed iterations before aborting fitting
 *  \param d_Parameters pointer to result array of fitted parameters, x, y coords, etc.
 *  \param d_CRLBs pointer to result array of Cramer-Rao lower bound estimates 
 *  \param d_LogLikelihood pointer to result array of loglikelihood estimates of fitting
 *  \param Nfits number of subregions to fit
 */
	kernel_MLEFit_LM_sigmaxy_EMCCD<<<dimGrid, dimBlock>>>(d_data, PSFSigma, sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits);
}


extern void kernel_splineMLEFit_z_EMCCD_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data,const float *d_coeff, const int spline_xsize, const int spline_ysize, const int spline_zsize, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits) 
{
	kernel_splineMLEFit_z_EMCCD<<<dimGrid, dimBlock>>>(d_data, d_coeff, spline_xsize, spline_ysize, spline_zsize, sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits);
}






//
////sCMOS wrapper
////*******************************************************************************************
//extern void kernel_MLEFit_sCMOS_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma, const int sz, const int iterations, 
//        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits, const float *d_varim, const float *d_gainim) 
//{
///*!
// *  \brief Basic maximum likelihood estimator fit based kernel
// *  \param dimGrid number of blocks 
// *  \param dimBlock number of threads per block
// *  \param d_data an array of subregions to be processed copied into video memory
// *  \param PSFSigma the sigma value to use for the point spread function
// *  \param sz nxn size of the subregion
// *  \param iterations maximum allowed iterations before aborting fitting
// *  \param d_Parameters pointer to result array of fitted parameters, x, y coords, etc.
// *  \param d_CRLBs pointer to result array of Cramer-Rao lower bound estimates 
// *  \param d_LogLikelihood pointer to result array of loglikelihood estimates of fitting
// *  \param Nfits number of subregions to fit
// */
//
//	//kernel_MLEFit_LM_sCMOS<<<dimGrid, dimBlock>>>(d_data, PSFSigma, sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits,d_varim, d_gainim);
//}
//
////*******************************************************************************************
//extern void kernel_MLEFit_sigma_sCMOS_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma, const int sz, const int iterations, 
//        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits, const float *d_varim, const float *d_gainim) 
//{
///*!
// *  \brief Basic maximum likelihood estimator fit based kernel
// *  \param dimGrid number of blocks 
// *  \param dimBlock number of threads per block
// *  \param d_data an array of subregions to be processed copied into video memory
// *  \param PSFSigma the sigma value to use for the point spread function
// *  \param sz nxn size of the subregion
// *  \param iterations maximum allowed iterations before aborting fitting
// *  \param d_Parameters pointer to result array of fitted parameters, x, y coords, etc.
// *  \param d_CRLBs pointer to result array of Cramer-Rao lower bound estimates 
// *  \param d_LogLikelihood pointer to result array of loglikelihood estimates of fitting
// *  \param Nfits number of subregions to fit
// */
//	//kernel_MLEFit_LM_sCMOS<<<dimGrid, dimBlock>>>(d_data, PSFSigma, sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits,d_varim, d_gainim);
//}
//
////*******************************************************************************************
//extern void kernel_MLEFit_z_sCMOS_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma_x, const float Ax, const float Ay, const float Bx, 
//		const float By, const float gamma, const float d, const float PSFSigma_y, const int sz, const int iterations, 
//        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits, const float *d_varim, const float *d_gainim) 
//{
///*!
// *  \brief Basic maximum likelihood estimator fit based kernel
// *  \param dimGrid number of blocks 
// *  \param dimBlock number of threads per block
// *  \param d_data an array of subregions to be processed copied into video memory
// *  \param PSFSigma_x the sigma value to use for the point spread function on the x axis
// *  \param Ax ???
// *  \param Ay ???
// *  \param Bx ???
// *  \param By ???
// *  \param gamma ???
// *  \param d ???
// *  \param PSFSigma_y the sigma value to use for the point spread function on the y axis
// *  \param sz nxn size of the subregion
// *  \param iterations maximum allowed iterations before aborting fitting
// *  \param d_Parameters pointer to result array of fitted parameters, x, y coords, etc.
// *  \param d_CRLBs pointer to result array of Cramer-Rao lower bound estimates 
// *  \param d_LogLikelihood pointer to result array of loglikelihood estimates of fitting
// *  \param Nfits number of subregions to fit
// */
//	//kernel_MLEFit_LM_sCMOS<<<dimGrid, dimBlock>>>(d_data, PSFSigma_x, Ax, Ay, Bx, By, gamma, d, PSFSigma_y, sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits,d_varim, d_gainim);
//}
//
////*******************************************************************************************
//extern void kernel_MLEFit_sigmaxy_sCMOS_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma, const int sz, const int iterations, 
//        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits, const float *d_varim, const float *d_gainim) 
//{
//	/*!
// *  \brief Basic maximum likelihood estimator fit based kernel
// *  \param dimGrid number of blocks 
// *  \param dimBlock number of threads per block
// *  \param d_data an array of subregions to be processed copied into video memory
// *  \param PSFSigma the sigma value to use for the point spread function
// *  \param sz nxn size of the subregion
// *  \param iterations maximum allowed iterations before aborting fitting
// *  \param d_Parameters pointer to result array of fitted parameters, x, y coords, etc.
// *  \param d_CRLBs pointer to result array of Cramer-Rao lower bound estimates 
// *  \param d_LogLikelihood pointer to result array of loglikelihood estimates of fitting
// *  \param Nfits number of subregions to fit
// */
//	//kernel_MLEFit_LM_sCMOS<<<dimGrid, dimBlock>>>(d_data, PSFSigma, sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits,d_varim, d_gainim);
//}


//extern void kernel_splineMLEFit_z_sCMOS_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data,const float *d_coeff, const int spline_xsize, const int spline_ysize, const int spline_zsize, const int sz, const int iterations, 
//        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits, const float *d_varim, const float *d_gainim) 
//{
//	kernel_splineMLEFit_z_sCMOS<<<dimGrid, dimBlock>>>(d_data, d_coeff, spline_xsize, spline_ysize, spline_zsize, sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits, d_varim, d_gainim);
//}