/*!
 * \file wrapper.cu
 * \author Keith Lidke
 * \date January 10, 2010
 * \brief Wrap the Cuda kernel calls as standard external C functions.  This allows the kernels to be
 * called without doing anything special in the C code and simplifies building the code.
 */

#include "GPUgaussMLEv2.h"

//*******************************************************************************************
extern "C" void kernel_MLEFit_wrapper(dim3 dimGrid, dim3 dimBlock, float *d_data, float PSFSigma, int sz, int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,int Nfits) 
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

	kernel_MLEFit<<<dimGrid, dimBlock>>>(d_data, PSFSigma, sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits);
}

//*******************************************************************************************
extern "C" void kernel_MLEFit_sigma_wrapper(dim3 dimGrid, dim3 dimBlock, float *d_data, float PSFSigma, int sz, int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,int Nfits) 
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
	kernel_MLEFit_sigma<<<dimGrid, dimBlock>>>(d_data, PSFSigma, sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits);
}

//*******************************************************************************************
extern "C" void kernel_MLEFit_z_wrapper(dim3 dimGrid, dim3 dimBlock, float *d_data, float PSFSigma_x, float Ax, float Ay, float Bx, 
		float By, float gamma, float d, float PSFSigma_y, int sz, int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,int Nfits) 
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
	kernel_MLEFit_z<<<dimGrid, dimBlock>>>(d_data, PSFSigma_x, Ax, Ay, Bx, By, gamma, d, PSFSigma_y, sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits);
}

//*******************************************************************************************
extern "C" void kernel_MLEFit_sigmaxy_wrapper(dim3 dimGrid, dim3 dimBlock, float *d_data, float PSFSigma, int sz, int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,int Nfits) 
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
	kernel_MLEFit_sigmaxy<<<dimGrid, dimBlock>>>(d_data, PSFSigma, sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits);
}