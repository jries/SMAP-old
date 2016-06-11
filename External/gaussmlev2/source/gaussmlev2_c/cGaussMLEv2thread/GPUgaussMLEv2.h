/*!
 * \file GPUgaussMLEv2.h
 * \author Keith Lidke
 * \date January 10, 2010
 * \brief Prototypes for the actual Cuda kernels.  
 */

#ifndef GPUGAUSSMLEV2_H
#define GPUGAUSSMLEV2_H

//void kernel_MLEFit_noshared(const int subregion, const float *d_data, const float PSFSigma, const int sz, const int iterations, 
//        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits);

void kernel_MLEFit(const int subregion, const float *d_data, const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits);

void kernel_MLEFit_sigma(const int subregion, const float *d_data, const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits);

void kernel_MLEFit_z(const int subregion, const float *d_data, const float PSFSigma_x, const float Ax, const float Ay, const float Bx, 
		const float By, const float gamma, const float d, const float PSFSigma_y, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits);

void kernel_MLEFit_sigmaxy(const int subregion, const float *d_data, const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits);

#endif