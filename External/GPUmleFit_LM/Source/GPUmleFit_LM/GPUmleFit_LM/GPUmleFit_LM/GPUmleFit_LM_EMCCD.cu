#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"
#include "cuda_runtime.h"
#include "definitions.h"
#include "MatInvLib.h"

#include "GPUsplineLib.cuh"
#include "GPUgaussLib.cuh"
#include "GPUmleFit_LM_EMCCD.h"


//#include "GPUsplineMLE.h"

__global__ void kernel_MLEFit_LM_EMCCD(const float *d_data,const float PSFSigma, const int sz, const int iterations, 
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits){

		const int NV=NV_P;
		float M[NV*NV],Diag[NV], Minv[NV*NV];
		int tx = threadIdx.x;
		int bx = blockIdx.x;
		int BlockSize = blockDim.x;
		int ii, jj, kk, ll, l, m, i;


		float model, data;
		float Div;
		float PSFy, PSFx;

		float newTheta[NV],oldTheta[NV];
		float newLambda = 1.0, oldLambda = 1.0;
		float newSign[NV] = {0}, oldSign[NV] = {0};
		float newUpdate[NV] = {0},oldUpdate[NV] = {0};
		float newClamp[NV]={1.0,1.0,100,20}, oldClamp[NV]={1.0,1.0,100,20};
		float newDudt[NV] ={0};

		float newErr = 1e13, oldErr = 1e13;

		float jacobian[NV]={0};
		float hessian[NV*NV]={0};
		float t1,t2;

		float Nmax;
		float temp;
		int info;
		float L[NV*NV] = {0}, U[NV*NV] = {0};


		//Prevent read/write past end of array
		if ((bx*BlockSize+tx)>=Nfits) return;

		for (ii=0;ii<NV*NV;ii++)M[ii]=0;
		for (ii=0;ii<NV*NV;ii++)Minv[ii]=0;

		//copy in data
		const float *s_data = d_data+(sz*sz*bx*BlockSize+sz*sz*tx);

		//initial values
		kernel_CenterofMass2D(sz, s_data, &newTheta[0], &newTheta[1]);
		kernel_GaussFMaxMin2D(sz, PSFSigma, s_data, &Nmax, &newTheta[3]);
		newTheta[2]=max(0.0, (Nmax-newTheta[3])*2*pi*PSFSigma*PSFSigma);


		for (ii=0;ii<NV;ii++)oldTheta[ii]=newTheta[ii];

		//updateFitValues3D
		newErr = 0;
		memset(jacobian,0,NV*sizeof(float));
		memset(hessian,0,NV*NV*sizeof(float));
		for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
			kernel_DerivativeGauss2D(ii,jj,PSFSigma,newTheta,newDudt,&model);
			data=s_data[sz*jj+ii];

			if (data>0)
				newErr = newErr + 2*((model-data)-data*log(model/data));
			else if (data==0)
				newErr = newErr + 2*model;

			t1 = 1-data/model;
			for (l=0;l<NV;l++){
				jacobian[l]+=t1*newDudt[l];
			}

			t2 = data/pow(model,2);
			for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
				hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
				hessian[m*NV+l] = hessian[l*NV+m];
			}
		}
		//addPeak

		//copyFitData

		for (kk=0;kk<iterations;kk++) {//main iterative loop

			if(abs((newErr-oldErr)/newErr)<TOLERANCE){
				//newStatus = CONVERGED;
				break;
			}
			else{
				if(newErr>1.5*oldErr){
					//copy Fitdata

					for (i=0;i<NV;i++){
						newSign[i]=oldSign[i];
						newClamp[i]=oldClamp[i];
						newTheta[i]=oldTheta[i];
					}
					newLambda = oldLambda;
					newErr = oldErr;

					newLambda = 10*newLambda;
				}
				else if(newErr<oldErr){
					if (newLambda>1){
						newLambda = newLambda*0.8;
					}
					else if(newLambda<1){
						newLambda = 1;
					}
				}


				for (i=0;i<NV;i++){
					hessian[i*NV+i]=hessian[i*NV+i]*newLambda;
				}
				memset(L,0,NV*sizeof(float));
				memset(U,0,NV*sizeof(float));
				info = kernel_cholesky(hessian,NV,L,U);
				if (info ==0){
					kernel_luEvaluate(L,U,jacobian,NV,newUpdate);
					//copyFitData
					for (i=0;i<NV;i++){
						oldSign[i]=newSign[i];
						oldClamp[i]=newClamp[i];

						oldTheta[i]=newTheta[i];
					}
					oldLambda = newLambda;
					oldErr=newErr;


					//updatePeakParameters
					for (ll=0;ll<NV;ll++){
						if (newSign[ll]!=0){
							if (newSign[ll]==1&&newUpdate[ll]<0){
								newClamp[ll]=newClamp[ll]*0.5;
							}
							else if (newSign[ll]==-1&&newUpdate[ll]>0){
								newClamp[ll] = newClamp[ll]*0.5;
							}
						}

						if (newUpdate[ll]>0){
							newSign[ll]=1;
						}
						else{
							newSign[ll]=-1;
						}

						newTheta[ll] = newTheta[ll]-newUpdate[ll]/(1+abs(newUpdate[ll]/newClamp[ll]));
					}

					newTheta[0] = max(newTheta[0],(float(sz)-1)/2-2);
					newTheta[0] = min(newTheta[0],(float(sz)-1)/2+2);
					newTheta[1] = max(newTheta[1],(float(sz)-1)/2-2);
					newTheta[1] = min(newTheta[1],(float(sz)-1)/2+2);
					newTheta[2] = max(newTheta[2],1.0);
					newTheta[3] = max(newTheta[3],0.01);
					


					newErr = 0;
					memset(jacobian,0,NV*sizeof(float));
					memset(hessian,0,NV*NV*sizeof(float));
					for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
						//calculating derivatives
						kernel_DerivativeGauss2D(ii,jj,PSFSigma,newTheta,newDudt,&model);
						data=s_data[sz*jj+ii];			

						if (data>0)
							newErr = newErr + 2*((model-data)-data*log(model/data));
						else if (data==0)
							newErr = newErr + 2*model;

						t1 = 1-data/model;
						for (l=0;l<NV;l++){
							jacobian[l]+=t1*newDudt[l];
						}

						t2 = data/pow(model,2);
						for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
							hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
							hessian[m*NV+l] = hessian[l*NV+m];
						}
					}
				}
				else
				{
					newLambda = 10*newLambda;
				}
			}
		}

		// Calculating the CRLB and LogLikelihood
		Div=0.0f;
		for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
			//need to check why don't use newTheta[4] instead of PSFSigma!!!
			kernel_DerivativeGauss2D(ii,jj,PSFSigma,newTheta,newDudt,&model);
			data=s_data[sz*jj+ii];

			//Building the Fisher Information Matrix
			for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
				M[kk*NV+ll]+= newDudt[ll]*newDudt[kk]/model;
				M[ll*NV+kk]=M[kk*NV+ll];
			}

			//LogLikelyhood
			if (model>0)
				if (data>0)Div+=data*log(model)-model-data*log(data)+data;
				else
					Div+=-model;
		}

		// Matrix inverse (CRLB=F^-1) and output assigments
		kernel_MatInvN(M, Minv, Diag, NV);
		//write to global arrays
		for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+BlockSize*bx+tx]=newTheta[kk];
		for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+BlockSize*bx+tx]=Diag[kk];
		d_LogLikelihood[BlockSize*bx+tx] = Div;

		return;
}

//*********************************************************************************************************************************************
//Continue to optimize Hazen's workflow  change xc,xStart,yc,yStart calculation. Remove some redundency in kernel_computeDelta3D and kernal_fAt3D

__global__ void kernel_MLEFit_LM_Sigma_EMCCD(const float *d_data,const float PSFSigma, const int sz, const int iterations, 
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits){

		const int NV=NV_PS;
		float M[NV*NV],Diag[NV], Minv[NV*NV];
		int tx = threadIdx.x;
		int bx = blockIdx.x;
		int BlockSize = blockDim.x;
		int ii, jj, kk, ll, l, m, i;
		//int xstart, ystart, zstart, xi, yi;


		float model, data;
		float Div;
		float PSFy, PSFx;

		float newTheta[NV],oldTheta[NV];
		float newLambda = 1.0, oldLambda = 1.0;
		float newSign[NV] = {0}, oldSign[NV] = {0};
		float newUpdate[NV] = {0},oldUpdate[NV] = {0};
		float newClamp[NV]={1.0,1.0,100,20,0.5}, oldClamp[NV]={1.0,1.0,100,20,0.5};
		float newDudt[NV] ={0};

		float newErr = 1e13, oldErr = 1e13;

		float jacobian[NV]={0};
		float hessian[NV*NV]={0};
		float t1,t2;

		float Nmax;
		float temp;
		int info;
		float L[NV*NV] = {0}, U[NV*NV] = {0};


		//Prevent read/write past end of array
		if ((bx*BlockSize+tx)>=Nfits) return;

		for (ii=0;ii<NV*NV;ii++)M[ii]=0;
		for (ii=0;ii<NV*NV;ii++)Minv[ii]=0;

		//copy in data
		const float *s_data = d_data+(sz*sz*bx*BlockSize+sz*sz*tx);

		//initial values
		kernel_CenterofMass2D(sz, s_data, &newTheta[0], &newTheta[1]);
		kernel_GaussFMaxMin2D(sz, PSFSigma, s_data, &Nmax, &newTheta[3]);
		newTheta[2]=max(0.0, (Nmax-newTheta[3])*2*pi*PSFSigma*PSFSigma);
		newTheta[4]=PSFSigma;

		for (ii=0;ii<NV;ii++)oldTheta[ii]=newTheta[ii];

		//updateFitValues3D
		newErr = 0;
		memset(jacobian,0,NV*sizeof(float));
		memset(hessian,0,NV*NV*sizeof(float));
		for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
			kernel_DerivativeGauss2D_sigma(ii,jj,newTheta,newDudt,&model);
			data=s_data[sz*jj+ii];

			if (data>0)
				newErr = newErr + 2*((model-data)-data*log(model/data));
			else if (data==0)
				newErr = newErr + 2*model;

			t1 = 1-data/model;
			for (l=0;l<NV;l++){
				jacobian[l]+=t1*newDudt[l];
			}

			t2 = data/pow(model,2);
			for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
				hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
				hessian[m*NV+l] = hessian[l*NV+m];
			}
		}
		//addPeak

		//copyFitData

		for (kk=0;kk<iterations;kk++) {//main iterative loop

			if(abs((newErr-oldErr)/newErr)<TOLERANCE){
				//newStatus = CONVERGED;
				break;
			}
			else{
				if(newErr>1.5*oldErr){
					//copy Fitdata

					for (i=0;i<NV;i++){
						newSign[i]=oldSign[i];
						newClamp[i]=oldClamp[i];
						newTheta[i]=oldTheta[i];
					}
					newLambda = oldLambda;
					newErr = oldErr;

					newLambda = 10*newLambda;
				}
				else if(newErr<oldErr){
					if (newLambda>1){
						newLambda = newLambda*0.8;
					}
					else if(newLambda<1){
						newLambda = 1;
					}
				}


				for (i=0;i<NV;i++){
					hessian[i*NV+i]=hessian[i*NV+i]*newLambda;
				}
				memset(L,0,NV*sizeof(float));
				memset(U,0,NV*sizeof(float));
				info = kernel_cholesky(hessian,NV,L,U);
				if (info ==0){
					kernel_luEvaluate(L,U,jacobian,NV,newUpdate);
					//copyFitData
					for (i=0;i<NV;i++){
						oldSign[i]=newSign[i];
						oldClamp[i]=newClamp[i];

						oldTheta[i]=newTheta[i];
					}
					oldLambda = newLambda;
					oldErr=newErr;


					//updatePeakParameters
					for (ll=0;ll<NV;ll++){
						if (newSign[ll]!=0){
							if (newSign[ll]==1&&newUpdate[ll]<0){
								newClamp[ll]=newClamp[ll]*0.5;
							}
							else if (newSign[ll]==-1&&newUpdate[ll]>0){
								newClamp[ll] = newClamp[ll]*0.5;
							}
						}

						if (newUpdate[ll]>0){
							newSign[ll]=1;
						}
						else{
							newSign[ll]=-1;
						}

						newTheta[ll] = newTheta[ll]-newUpdate[ll]/(1+abs(newUpdate[ll]/newClamp[ll]));
					}

					newTheta[0] = max(newTheta[0],(float(sz)-1)/2-2);
					newTheta[0] = min(newTheta[0],(float(sz)-1)/2+2);
					newTheta[1] = max(newTheta[1],(float(sz)-1)/2-2);
					newTheta[1] = min(newTheta[1],(float(sz)-1)/2+2);
					newTheta[2] = max(newTheta[2],1.0);
					newTheta[3] = max(newTheta[3],0.01);
					newTheta[4] = max(newTheta[4],0.0);
					newTheta[4] = min(newTheta[4],sz/2.0f);


					newErr = 0;
					memset(jacobian,0,NV*sizeof(float));
					memset(hessian,0,NV*NV*sizeof(float));
					for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
						//calculating derivatives
						kernel_DerivativeGauss2D_sigma(ii,jj,newTheta,newDudt,&model);
						data=s_data[sz*jj+ii];			

						if (data>0)
							newErr = newErr + 2*((model-data)-data*log(model/data));
						else if (data==0)
							newErr = newErr + 2*model;

						t1 = 1-data/model;
						for (l=0;l<NV;l++){
							jacobian[l]+=t1*newDudt[l];
						}

						t2 = data/pow(model,2);
						for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
							hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
							hessian[m*NV+l] = hessian[l*NV+m];
						}
					}
				}
				else
				{
					newLambda = 10*newLambda;
				}
			}
		}

		// Calculating the CRLB and LogLikelihood
		Div=0.0f;
		for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
			//need to check why don't use newTheta[4] instead of PSFSigma!!!
			//PSFx=kernel_IntGauss1D(ii, newTheta[0], PSFSigma);
			//PSFy=kernel_IntGauss1D(jj, newTheta[1], PSFSigma);

			//model=newTheta[3]+newTheta[2]*PSFx*PSFy;
			//data=s_data[sz*jj+ii];

			////calculating derivatives
			//kernel_DerivativeIntGauss1D(ii, newTheta[0], newTheta[4], newTheta[2], PSFy, &newDudt[0], NULL);
			//kernel_DerivativeIntGauss1D(jj, newTheta[1], newTheta[4], newTheta[2], PSFx, &newDudt[1], NULL);
			//kernel_DerivativeIntGauss2DSigma(ii, jj, newTheta[0], newTheta[1], newTheta[4], newTheta[2], PSFx, PSFy, &newDudt[4], NULL);
			//newDudt[2] = PSFx*PSFy;
			//newDudt[3] = 1.0f;
			kernel_DerivativeGauss2D_sigma(ii,jj,newTheta,newDudt,&model);
			data=s_data[sz*jj+ii];

			//Building the Fisher Information Matrix
			for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
				M[kk*NV+ll]+= newDudt[ll]*newDudt[kk]/model;
				M[ll*NV+kk]=M[kk*NV+ll];
			}

			//LogLikelyhood
			if (model>0)
				if (data>0)Div+=data*log(model)-model-data*log(data)+data;
				else
					Div+=-model;
		}

		// Matrix inverse (CRLB=F^-1) and output assigments
		kernel_MatInvN(M, Minv, Diag, NV);
		//write to global arrays
		for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+BlockSize*bx+tx]=newTheta[kk];
		for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+BlockSize*bx+tx]=Diag[kk];
		d_LogLikelihood[BlockSize*bx+tx] = Div;

		return;
}

//*********************************************************************************************************************************************
//Continue to optimize Hazen's workflow  change xc,xStart,yc,yStart calculation. Remove some redundency in kernel_computeDelta3D and kernal_fAt3D

__global__ void kernel_MLEFit_LM_z_EMCCD(const float *d_data, const float PSFSigma_x, const float Ax, const float Ay, const float Bx, 
	const float By, const float gamma, const float d, const float PSFSigma_y, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits){

		const int NV=NV_PZ;
		float M[NV*NV],Diag[NV], Minv[NV*NV];
		int tx = threadIdx.x;
		int bx = blockIdx.x;
		int BlockSize = blockDim.x;
		int ii, jj, kk, ll, l, m, i;
		//int xstart, ystart, zstart, xi, yi;


		float model, data;
		float Div;
		float PSFy, PSFx;

		float newTheta[NV],oldTheta[NV];
		float newLambda = 1.0, oldLambda = 1.0;
		float newSign[NV] = {0}, oldSign[NV] = {0};
		float newUpdate[NV] = {0},oldUpdate[NV] = {0};
		float newClamp[NV]={1.0,1.0,100,20,2}, oldClamp[NV]={1.0,1.0,100,20,2};
		float newDudt[NV] ={0};

		float newErr = 1e13, oldErr = 1e13;

		float jacobian[NV]={0};
		float hessian[NV*NV]={0};
		float t1,t2;

		float Nmax;
		float temp;
		int info;
		float L[NV*NV] = {0}, U[NV*NV] = {0};


		//Prevent read/write past end of array
		if ((bx*BlockSize+tx)>=Nfits) return;

		for (ii=0;ii<NV*NV;ii++)M[ii]=0;
		for (ii=0;ii<NV*NV;ii++)Minv[ii]=0;

		//copy in data
		const float *s_data = d_data+(sz*sz*bx*BlockSize+sz*sz*tx);

		//initial values
		kernel_CenterofMass2D(sz, s_data, &newTheta[0], &newTheta[1]);
		kernel_GaussFMaxMin2D(sz, PSFSigma_x, s_data, &Nmax, &newTheta[3]);
		newTheta[2]=max(0.0, (Nmax-newTheta[3])*2*pi*PSFSigma_x*PSFSigma_y*sqrt(2.0f));
		newTheta[4]=0;

		for (ii=0;ii<NV;ii++)oldTheta[ii]=newTheta[ii];

		//updateFitValues3D
		newErr = 0;
		memset(jacobian,0,NV*sizeof(float));
		memset(hessian,0,NV*NV*sizeof(float));
		for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
			 kernel_DerivativeIntGauss2Dz(ii, jj, newTheta, PSFSigma_x,PSFSigma_y, Ax,Ay,Bx,By, gamma, d, &PSFx, &PSFy, newDudt, NULL,&model);
			data=s_data[sz*jj+ii];

			if (data>0)
				newErr = newErr + 2*((model-data)-data*log(model/data));
			else if (data==0)
				newErr = newErr + 2*model;

			t1 = 1-data/model;
			for (l=0;l<NV;l++){
				jacobian[l]+=t1*newDudt[l];
			}

			t2 = data/pow(model,2);
			for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
				hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
				hessian[m*NV+l] = hessian[l*NV+m];
			}
		}
		//addPeak

		//copyFitData

		for (kk=0;kk<iterations;kk++) {//main iterative loop

			if(abs((newErr-oldErr)/newErr)<TOLERANCE){
				//newStatus = CONVERGED;
				break;
			}
			else{
				if(newErr>1.5*oldErr){
					//copy Fitdata

					for (i=0;i<NV;i++){
						newSign[i]=oldSign[i];
						newClamp[i]=oldClamp[i];
						newTheta[i]=oldTheta[i];
					}
					newLambda = oldLambda;
					newErr = oldErr;

					newLambda = 10*newLambda;
				}
				else if(newErr<oldErr){
					if (newLambda>1){
						newLambda = newLambda*0.8;
					}
					else if(newLambda<1){
						newLambda = 1;
					}
				}


				for (i=0;i<NV;i++){
					hessian[i*NV+i]=hessian[i*NV+i]*newLambda;
				}
				memset(L,0,NV*sizeof(float));
				memset(U,0,NV*sizeof(float));
				info = kernel_cholesky(hessian,NV,L,U);
				if (info ==0){
					kernel_luEvaluate(L,U,jacobian,NV,newUpdate);
					//copyFitData
					for (i=0;i<NV;i++){
						oldSign[i]=newSign[i];
						oldClamp[i]=newClamp[i];

						oldTheta[i]=newTheta[i];
					}
					oldLambda = newLambda;
					oldErr=newErr;


					//updatePeakParameters
					for (ll=0;ll<NV;ll++){
						if (newSign[ll]!=0){
							if (newSign[ll]==1&&newUpdate[ll]<0){
								newClamp[ll]=newClamp[ll]*0.5;
							}
							else if (newSign[ll]==-1&&newUpdate[ll]>0){
								newClamp[ll] = newClamp[ll]*0.5;
							}
						}

						if (newUpdate[ll]>0){
							newSign[ll]=1;
						}
						else{
							newSign[ll]=-1;
						}

						newTheta[ll] = newTheta[ll]-newUpdate[ll]/(1+abs(newUpdate[ll]/newClamp[ll]));
					}

					newTheta[0] = max(newTheta[0],(float(sz)-1)/2-2);
					newTheta[0] = min(newTheta[0],(float(sz)-1)/2+2);
					newTheta[1] = max(newTheta[1],(float(sz)-1)/2-2);
					newTheta[1] = min(newTheta[1],(float(sz)-1)/2+2);
					newTheta[2] = max(newTheta[2],1.0);
					newTheta[3] = max(newTheta[3],0.01);
					//newTheta[4] = max(newTheta[4],0.0);
					//newTheta[4] = min(newTheta[4],sz/2.0f);


					newErr = 0;
					memset(jacobian,0,NV*sizeof(float));
					memset(hessian,0,NV*NV*sizeof(float));
					for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
						//calculating derivatives
						kernel_DerivativeGauss2D_sigma(ii,jj,newTheta,newDudt,&model);
						data=s_data[sz*jj+ii];			

						if (data>0)
							newErr = newErr + 2*((model-data)-data*log(model/data));
						else if (data==0)
							newErr = newErr + 2*model;

						t1 = 1-data/model;
						for (l=0;l<NV;l++){
							jacobian[l]+=t1*newDudt[l];
						}

						t2 = data/pow(model,2);
						for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
							hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
							hessian[m*NV+l] = hessian[l*NV+m];
						}
					}
				}
				else
				{
					newLambda = 10*newLambda;
				}
			}
		}

		// Calculating the CRLB and LogLikelihood
		Div=0.0f;
		for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
			//need to check why don't use newTheta[4] instead of PSFSigma!!!
		    kernel_DerivativeIntGauss2Dz(ii, jj, newTheta, PSFSigma_x,PSFSigma_y, Ax,Ay,Bx,By, gamma, d, &PSFx, &PSFy, newDudt, NULL,&model);
			data=s_data[sz*jj+ii];

			//Building the Fisher Information Matrix
			for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
				M[kk*NV+ll]+= newDudt[ll]*newDudt[kk]/model;
				M[ll*NV+kk]=M[kk*NV+ll];
			}

			//LogLikelyhood
			if (model>0)
				if (data>0)Div+=data*log(model)-model-data*log(data)+data;
				else
					Div+=-model;
		}

		// Matrix inverse (CRLB=F^-1) and output assigments
		kernel_MatInvN(M, Minv, Diag, NV);
		//write to global arrays
		for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+BlockSize*bx+tx]=newTheta[kk];
		for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+BlockSize*bx+tx]=Diag[kk];
		d_LogLikelihood[BlockSize*bx+tx] = Div;

		return;
}

//*********************************************************************************************************************************************
//Continue to optimize Hazen's workflow  change xc,xStart,yc,yStart calculation. Remove some redundency in kernel_computeDelta3D and kernal_fAt3D

__global__ void kernel_MLEFit_LM_sigmaxy_EMCCD(const float *d_data, const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits){

		const int NV=NV_PS2;
		float M[NV*NV],Diag[NV], Minv[NV*NV];
		int tx = threadIdx.x;
		int bx = blockIdx.x;
		int BlockSize = blockDim.x;
		int ii, jj, kk, ll, l, m, i;
		//int xstart, ystart, zstart, xi, yi;


		float model, data;
		float Div;
		float PSFy, PSFx;

		float newTheta[NV],oldTheta[NV];
		float newLambda = 1.0, oldLambda = 1.0;
		float newSign[NV] = {0}, oldSign[NV] = {0};
		float newUpdate[NV] = {0},oldUpdate[NV] = {0};
		float newClamp[NV]={1.0,1.0,100,20,0.5,0.5}, oldClamp[NV]={1.0,1.0,100,20,0.5,0.5};
		float newDudt[NV] ={0};

		float newErr = 1e13, oldErr = 1e13;

		float jacobian[NV]={0};
		float hessian[NV*NV]={0};
		float t1,t2;

		float Nmax;
		float temp;
		int info;
		float L[NV*NV] = {0}, U[NV*NV] = {0};


		//Prevent read/write past end of array
		if ((bx*BlockSize+tx)>=Nfits) return;

		for (ii=0;ii<NV*NV;ii++)M[ii]=0;
		for (ii=0;ii<NV*NV;ii++)Minv[ii]=0;

		//copy in data
		const float *s_data = d_data+(sz*sz*bx*BlockSize+sz*sz*tx);

		//initial values
		kernel_CenterofMass2D(sz, s_data, &newTheta[0], &newTheta[1]);
		kernel_GaussFMaxMin2D(sz, PSFSigma, s_data, &Nmax, &newTheta[3]);
		newTheta[2]=max(0.0, (Nmax-newTheta[3])*2*pi*PSFSigma*PSFSigma);
		newTheta[4]=PSFSigma;
		newTheta[5]=PSFSigma;

		for (ii=0;ii<NV;ii++)oldTheta[ii]=newTheta[ii];

		//updateFitValues3D
		newErr = 0;
		memset(jacobian,0,NV*sizeof(float));
		memset(hessian,0,NV*NV*sizeof(float));
		for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
			kernel_DerivativeGauss2D_sigmaxy( ii,  jj, PSFSigma, newTheta, newDudt, &model);
			data=s_data[sz*jj+ii];

			if (data>0)
				newErr = newErr + 2*((model-data)-data*log(model/data));
			else if (data==0)
				newErr = newErr + 2*model;

			t1 = 1-data/model;
			for (l=0;l<NV;l++){
				jacobian[l]+=t1*newDudt[l];
			}

			t2 = data/pow(model,2);
			for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
				hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
				hessian[m*NV+l] = hessian[l*NV+m];
			}
		}
		//addPeak

		//copyFitData

		for (kk=0;kk<iterations;kk++) {//main iterative loop

			if(abs((newErr-oldErr)/newErr)<TOLERANCE){
				//newStatus = CONVERGED;
				break;
			}
			else{
				if(newErr>1.5*oldErr){
					//copy Fitdata

					for (i=0;i<NV;i++){
						newSign[i]=oldSign[i];
						newClamp[i]=oldClamp[i];
						newTheta[i]=oldTheta[i];
					}
					newLambda = oldLambda;
					newErr = oldErr;

					newLambda = 10*newLambda;
				}
				else if(newErr<oldErr){
					if (newLambda>1){
						newLambda = newLambda*0.8;
					}
					else if(newLambda<1){
						newLambda = 1;
					}
				}


				for (i=0;i<NV;i++){
					hessian[i*NV+i]=hessian[i*NV+i]*newLambda;
				}
				memset(L,0,NV*sizeof(float));
				memset(U,0,NV*sizeof(float));
				info = kernel_cholesky(hessian,NV,L,U);
				if (info ==0){
					kernel_luEvaluate(L,U,jacobian,NV,newUpdate);
					//copyFitData
					for (i=0;i<NV;i++){
						oldSign[i]=newSign[i];
						oldClamp[i]=newClamp[i];

						oldTheta[i]=newTheta[i];
					}
					oldLambda = newLambda;
					oldErr=newErr;


					//updatePeakParameters
					for (ll=0;ll<NV;ll++){
						if (newSign[ll]!=0){
							if (newSign[ll]==1&&newUpdate[ll]<0){
								newClamp[ll]=newClamp[ll]*0.5;
							}
							else if (newSign[ll]==-1&&newUpdate[ll]>0){
								newClamp[ll] = newClamp[ll]*0.5;
							}
						}

						if (newUpdate[ll]>0){
							newSign[ll]=1;
						}
						else{
							newSign[ll]=-1;
						}

						newTheta[ll] = newTheta[ll]-newUpdate[ll]/(1+abs(newUpdate[ll]/newClamp[ll]));
					}

					newTheta[0] = max(newTheta[0],(float(sz)-1)/2-2);
					newTheta[0] = min(newTheta[0],(float(sz)-1)/2+2);
					newTheta[1] = max(newTheta[1],(float(sz)-1)/2-2);
					newTheta[1] = min(newTheta[1],(float(sz)-1)/2+2);
					newTheta[2] = max(newTheta[2],1.0);
					newTheta[3] = max(newTheta[3],0.01);
					newTheta[4] = max(newTheta[4],PSFSigma/10.0f);
					newTheta[5] = max(newTheta[5],PSFSigma/10.0f);


					newErr = 0;
					memset(jacobian,0,NV*sizeof(float));
					memset(hessian,0,NV*NV*sizeof(float));
					for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
						//calculating derivatives
						kernel_DerivativeGauss2D_sigma(ii,jj,newTheta,newDudt,&model);
						data=s_data[sz*jj+ii];			

						if (data>0)
							newErr = newErr + 2*((model-data)-data*log(model/data));
						else if (data==0)
							newErr = newErr + 2*model;

						t1 = 1-data/model;
						for (l=0;l<NV;l++){
							jacobian[l]+=t1*newDudt[l];
						}

						t2 = data/pow(model,2);
						for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
							hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
							hessian[m*NV+l] = hessian[l*NV+m];
						}
					}
				}
				else
				{
					newLambda = 10*newLambda;
				}
			}
		}

		// Calculating the CRLB and LogLikelihood
		Div=0.0f;
		for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
			//need to check why don't use newTheta[4] instead of PSFSigma!!!
			kernel_DerivativeGauss2D_sigmaxy( ii,  jj, PSFSigma, newTheta, newDudt, &model);
			data=s_data[sz*jj+ii];

			//Building the Fisher Information Matrix
			for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
				M[kk*NV+ll]+= newDudt[ll]*newDudt[kk]/model;
				M[ll*NV+kk]=M[kk*NV+ll];
			}

			//LogLikelyhood
			if (model>0)
				if (data>0)Div+=data*log(model)-model-data*log(data)+data;
				else
					Div+=-model;
		}

		// Matrix inverse (CRLB=F^-1) and output assigments
		kernel_MatInvN(M, Minv, Diag, NV);
		//write to global arrays
		for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+BlockSize*bx+tx]=newTheta[kk];
		for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+BlockSize*bx+tx]=Diag[kk];
		d_LogLikelihood[BlockSize*bx+tx] = Div;

		return;
}

//******************************************************************************************************
//Continue to optimize Hazen's workflow  change xc,xStart,yc,yStart calculation. Remove some redundency in kernel_computeDelta3D and kernal_fAt3D

__global__ void kernel_splineMLEFit_z_EMCCD(const float *d_data,const float *d_coeff, const int spline_xsize, const int spline_ysize, const int spline_zsize, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits){
	
   const int NV=NV_PSP;
    float M[NV*NV],Diag[NV], Minv[NV*NV];
    int tx = threadIdx.x;
    int bx = blockIdx.x;
    int BlockSize = blockDim.x;
    int ii, jj, kk, ll, l, m, i;
	int xstart, ystart, zstart, xi, yi;

	const float *s_coeff;
	s_coeff = d_coeff;

    float model, data;
    float Div;
    //float dudt[NV_PS];
    float newTheta[NV],oldTheta[NV];
	float newLambda = 1.0, oldLambda = 1.0;
	float newSign[NV] = {0}, oldSign[NV] = {0};
	float newUpdate[NV] = {0},oldUpdate[NV] = {0};
	float newClamp[NV]={1.0,1.0,100,20,2}, oldClamp[NV]={1.0,1.0,100,20,2};
	float newDudt[NV] ={0};

	float newErr = 1e12, oldErr = 1e13;

	float off;
	float jacobian[NV]={0};
	float hessian[NV*NV]={0};
	float t1,t2;

	float Nmax;
	float xc,yc,zc;
	float temp;
	float delta_f[64]={0}, delta_dxf[64], delta_ddxf[64], delta_dyf[64], delta_ddyf[64], delta_dzf[64], delta_ddzf[64];
	int info;
	float L[NV*NV] = {0}, U[NV*NV] = {0};

    
    //Prevent read/write past end of array
    if ((bx*BlockSize+tx)>=Nfits) return;
    
    for (ii=0;ii<NV*NV;ii++)M[ii]=0;
    for (ii=0;ii<NV*NV;ii++)Minv[ii]=0;

    //copy in data
      const float *s_data = d_data+(sz*sz*bx*BlockSize+sz*sz*tx);
	  //const float *s_varim = d_varim+(sz*sz*bx*BlockSize+sz*sz*tx);
	//const float *s_gainim = d_gainim+(sz*sz*bx*BlockSize+sz*sz*tx);
    
    //initial values
    kernel_CenterofMass2D(sz, s_data, &newTheta[0], &newTheta[1]);
    kernel_GaussFMaxMin2D(sz, 1.5, s_data, &Nmax, &newTheta[3]);
    newTheta[2]=max(0.0, (Nmax-newTheta[3])*2*pi*1.5*1.5);
    newTheta[4]=float(spline_zsize)/2;

	for (ii=0;ii<NV;ii++)oldTheta[ii]=newTheta[ii];

	//updateFitValues3D
	xc = -2.0*((newTheta[0]-float(sz)/2)+0.5);
	yc = -2.0*((newTheta[1]-float(sz)/2)+0.5);

	off = (float(spline_xsize)+1.0-2*float(sz))/2;

	xstart = floor(xc);
	xc = xc-xstart;

	ystart = floor(yc);
	yc = yc-ystart;

	//zstart = floor(newTheta[4]);
	zstart = floor(newTheta[4]);
	zc = newTheta[4] -zstart;

	newErr = 0;
	memset(jacobian,0,NV*sizeof(float));
	memset(hessian,0,NV*NV*sizeof(float));
	kernel_computeDelta3D(xc, yc, zc, delta_f, delta_dxf, delta_dyf, delta_dzf);
	for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
		temp = kernal_fAt3D(2*ii+xstart+off,2*jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,s_coeff);
		model = newTheta[3]+newTheta[2]*temp;
		data=s_data[sz*jj+ii];
		//calculating derivatives

		newDudt[0] = -1*newTheta[2]*kernal_fAt3D(2*ii+xstart+off,2*jj+ystart+off,zstart, spline_xsize,spline_ysize,spline_zsize,delta_dxf,s_coeff);
		newDudt[1] = -1*newTheta[2]*kernal_fAt3D(2*ii+xstart+off,2*jj+ystart+off,zstart, spline_xsize,spline_ysize,spline_zsize,delta_dyf,s_coeff);
		newDudt[4] = newTheta[2]*kernal_fAt3D(2*ii+xstart+off,2*jj+ystart+off,zstart, spline_xsize,spline_ysize,spline_zsize,delta_dzf,s_coeff);
		newDudt[2] = temp;
		newDudt[3] = 1;

		if (data>0)
			newErr = newErr + 2*((model-data)-data*log(model/data));
		else if (data==0)
			newErr = newErr + 2*model;

		t1 = 1-data/model;
		for (l=0;l<NV;l++){
			jacobian[l]+=t1*newDudt[l];
		}

		t2 = data/pow(model,2);
		for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
			hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
			hessian[m*NV+l] = hessian[l*NV+m];
		}
	}
	//addPeak

	//copyFitData

	for (kk=0;kk<iterations;kk++) {//main iterative loop

			if(abs((newErr-oldErr)/newErr)<TOLERANCE){
				//newStatus = CONVERGED;
				break;
			}
			else{
				if(newErr>1.5*oldErr){
					//copy Fitdata
					
					for (i=0;i<NV;i++){
						newSign[i]=oldSign[i];
						newClamp[i]=oldClamp[i];
						newTheta[i]=oldTheta[i];
					}
					newLambda = oldLambda;
					newErr = oldErr;

					newLambda = 10*newLambda;
				}
				else if(newErr<oldErr){
					if (newLambda>1){
						newLambda = newLambda*0.8;
					}
					else if(newLambda<1){
						newLambda = 1;
					}
				}
				

				for (i=0;i<NV;i++){
					hessian[i*NV+i]=hessian[i*NV+i]*newLambda;
				}
				memset(L,0,NV*sizeof(float));
				memset(U,0,NV*sizeof(float));
				info = kernel_cholesky(hessian,NV,L,U);
				if (info ==0){
					kernel_luEvaluate(L,U,jacobian,NV,newUpdate);
                 //copyFitData
					for (i=0;i<NV;i++){
						oldSign[i]=newSign[i];
						oldClamp[i]=newClamp[i];

						oldTheta[i]=newTheta[i];
					}
					oldLambda = newLambda;
					oldErr=newErr;
					

					//updatePeakParameters
					for (ll=0;ll<NV;ll++){
						if (newSign[ll]!=0){
							if (newSign[ll]==1&&newUpdate[ll]<0){
								newClamp[ll]=newClamp[ll]*0.5;
							}
							else if (newSign[ll]==-1&&newUpdate[ll]>0){
								newClamp[ll] = newClamp[ll]*0.5;
							}
						}

						if (newUpdate[ll]>0){
							newSign[ll]=1;
						}
						else{
							newSign[ll]=-1;
						}

						newTheta[ll] = newTheta[ll]-newUpdate[ll]/(1+abs(newUpdate[ll]/newClamp[ll]));
					}

					newTheta[0] = max(newTheta[0],(float(sz)-1)/2-2);
					newTheta[0] = min(newTheta[0],(float(sz)-1)/2+2);
					newTheta[1] = max(newTheta[1],(float(sz)-1)/2-2);
					newTheta[1] = min(newTheta[1],(float(sz)-1)/2+2);
					newTheta[2] = max(newTheta[2],1.0);
					newTheta[3] = max(newTheta[3],0.01);
					newTheta[4] = max(newTheta[4],0.0);
					newTheta[4] = min(newTheta[4],float(spline_zsize));

					//updateFitValues3D
					xc = -2.0*((newTheta[0]-float(sz)/2)+0.5);
					yc = -2.0*((newTheta[1]-float(sz)/2)+0.5);

					off = (float(spline_xsize)+1.0-2*float(sz))/2;

					xstart = floor(xc);
					xc = xc-xstart;

					ystart = floor(yc);
					yc = yc-ystart;

					//zstart = floor(newTheta[4]);
					zstart = floor(newTheta[4]);
					zc = newTheta[4] -zstart;


					newErr = 0;
					memset(jacobian,0,NV*sizeof(float));
					memset(hessian,0,NV*NV*sizeof(float));
					kernel_computeDelta3D(xc, yc, zc, delta_f, delta_dxf, delta_dyf, delta_dzf);
					for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
						temp = kernal_fAt3D(2*ii+xstart+off,2*jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,s_coeff);
						model = newTheta[3]+newTheta[2]*temp;
						data=s_data[sz*jj+ii];
						//calculating derivatives

						newDudt[0] = -1*newTheta[2]*kernal_fAt3D(2*ii+xstart+off,2*jj+ystart+off,zstart, spline_xsize,spline_ysize,spline_zsize,delta_dxf,s_coeff);
						newDudt[1] = -1*newTheta[2]*kernal_fAt3D(2*ii+xstart+off,2*jj+ystart+off,zstart, spline_xsize,spline_ysize,spline_zsize,delta_dyf,s_coeff);
						newDudt[4] = newTheta[2]*kernal_fAt3D(2*ii+xstart+off,2*jj+ystart+off,zstart, spline_xsize,spline_ysize,spline_zsize,delta_dzf,s_coeff);
						newDudt[2] = temp;
						newDudt[3] = 1;

						if (data>0)
							newErr = newErr + 2*((model-data)-data*log(model/data));
						else if (data==0)
							newErr = newErr + 2*model;

						t1 = 1-data/model;
						for (l=0;l<NV;l++){
							jacobian[l]+=t1*newDudt[l];
						}

						t2 = data/pow(model,2);
						for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
							hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
							hessian[m*NV+l] = hessian[l*NV+m];
						}
					}
				}
				else
				{
					newLambda = 10*newLambda;
				}

					//copyFitdata	

			}


		
	}
    
    // Calculating the CRLB and LogLikelihood
	Div=0.0;

	xc = -2.0*((newTheta[0]-float(sz)/2)+0.5);
	yc = -2.0*((newTheta[1]-float(sz)/2)+0.5);

	off = (float(spline_xsize)+1.0-2*float(sz))/2;

	xstart = floor(xc);
	xc = xc-xstart;

	ystart = floor(yc);
	yc = yc-ystart;

	//zstart = floor(newTheta[4]);
	zstart = floor(newTheta[4]);
	zc = newTheta[4] -zstart;

	kernel_computeDelta3D(xc, yc, zc, delta_f, delta_dxf, delta_dyf, delta_dzf);

    for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
		temp = kernal_fAt3D(2*ii+xstart+off,2*jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,s_coeff);
		model = newTheta[3]+newTheta[2]*temp;
		data=s_data[sz*jj+ii];
		//calculating derivatives
		newDudt[0] = -1*newTheta[2]*kernal_fAt3D(2*ii+xstart+off,2*jj+ystart+off,zstart, spline_xsize,spline_ysize,spline_zsize,delta_dxf,s_coeff);
		newDudt[1] = -1*newTheta[2]*kernal_fAt3D(2*ii+xstart+off,2*jj+ystart+off,zstart, spline_xsize,spline_ysize,spline_zsize,delta_dyf,s_coeff);
		newDudt[4] = newTheta[2]*kernal_fAt3D(2*ii+xstart+off,2*jj+ystart+off,zstart, spline_xsize,spline_ysize,spline_zsize,delta_dzf,s_coeff);
		newDudt[2] = temp;
		newDudt[3] = 1;
        
        //Building the Fisher Information Matrix
        for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
            M[kk*NV+ll]+= newDudt[ll]*newDudt[kk]/model;
            M[ll*NV+kk]=M[kk*NV+ll];
        }
        
        //LogLikelyhood
        if (model>0)
            if (data>0)Div+=data*log(model)-model-data*log(data)+data;
            else
                Div+=-model;
    }
    
    // Matrix inverse (CRLB=F^-1) and output assigments
    kernel_MatInvN(M, Minv, Diag, NV);
  
    
    //write to global arrays
    for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+BlockSize*bx+tx]=newTheta[kk];
   for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+BlockSize*bx+tx]=Diag[kk];
   d_LogLikelihood[BlockSize*bx+tx] = Div;
	//d_LogLikelihood[BlockSize*bx+tx] = 1;
    
    
    return;
}



////*********************************************************************************************************************************************
////Continue to optimize Hazen's workflow  change xc,xStart,yc,yStart calculation. Remove some redundency in kernel_computeDelta3D and kernal_fAt3D
//
//__global__ void kernel_MLEFit_LM_EMCCD(const float *d_data,InitPara *para,float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,int fittype){
//
//	
//		int ii, jj, kk, ll, l, m, i;
//		
//		//const float PSFSigma= para->PSFSigma;
//		//const int sz = para->sz;
//		//const int iterations = para->iterations;
//		//const int Nfits = para->Nfits;
//	    //const int NV = para->NV;
//		const int NV = 5;
//		//float *newClamp=(float*)malloc(NV*sizeof(float));
//		//float *oldClamp=(float*)malloc(NV*sizeof(float));
//		//for (ii=0;ii<NV;ii++){
//		//	newClamp[ii]=para->clamp[ii];
//		//	oldClamp[ii]=para->clamp[ii];
//		//}
//		float newClamp[NV]={1.0,1.0,100,20,0.5}, oldClamp[NV]={1.0,1.0,100,20,0.5};
//		float PSFSigma = 1.5;
//		int sz = 7;
//		int iterations = 20;
//		int Nfits = 10000;
//
//		//const float PSFSigma_x=para->PSFSigma_x;
//		//const float PSFSigma_y=para->PSFSigma_y;
//		//const float Ax=para->Ax;
//		//const float Ay=para->Ay;
//		//const float Bx=para->Bx;
//		//const float By=para->By;
//		//const float gamma=para->gamma;
//		//const float d=para->d;
//
//		//float *d_coeff=para->coeff;
//		//int spline_xsize=para->spline_xsize;
//		//int spline_ysize=para->spline_ysize;
//		//int spline_zsize=para->spline_ysize;
//		//float newTheta[NV],oldTheta[NV];
//		float *newTheta=(float*)malloc(NV*sizeof(float));
//		float *oldTheta=(float*)malloc(NV*sizeof(float));
//
//
//		float Nmax;
//		int xstart, ystart, zstart, xi, yi;
//		float xc,yc,zc;
//		float off;
//		float delta_f[64]={0}, delta_dxf[64], delta_ddxf[64], delta_dyf[64], delta_ddyf[64], delta_dzf[64], delta_ddzf[64];
//		
//
//		int tx = threadIdx.x;
//		int bx = blockIdx.x;
//		int BlockSize = blockDim.x;
//		//Prevent read/write past end of array
//		if ((bx*BlockSize+tx)>=Nfits) return;
//		//copy in data
//		const float *s_data = d_data+(sz*sz*bx*BlockSize+sz*sz*tx);
//
//
//		//***************************************************************************************************************
//		// Parameter initilaization
//		switch (fittype){
//		case 1:
//			//float newClamp[NV]={1.0,1.0,100,20}, oldClamp[NV]={1.0,1.0,100,20};
//
//			//kernel_CenterofMass2D(sz, s_data, &newTheta[0], &newTheta[1]);
//			//kernel_GaussFMaxMin2D(sz, PSFSigma, s_data, &Nmax, &newTheta[3]);
//			//newTheta[2]=max(0.0f, (Nmax-newTheta[3])*2*pi*PSFSigma*PSFSigma);
//
//		case 2:
//			//float newClamp[NV]={1.0,1.0,100,20,0.5}, oldClamp[NV]={1.0,1.0,100,20,0.5};
//
//			kernel_CenterofMass2D(sz, s_data, &newTheta[0], &newTheta[1]);
//			kernel_GaussFMaxMin2D(sz, PSFSigma, s_data, &Nmax, &newTheta[3]);
//			newTheta[2]=max(0.0f, (Nmax-newTheta[3])*2*pi*PSFSigma*PSFSigma);
//			newTheta[4]=PSFSigma;
//
//		case 3:
//			//float newClamp[NV]={1.0,1.0,100,20,1.0}, oldClamp[NV]={1.0,1.0,100,20,1.0};
//
//			//kernel_CenterofMass2D(sz, s_data, &newTheta[0], &newTheta[1]);
//			//kernel_GaussFMaxMin2D(sz, PSFSigma_x, s_data, &Nmax, &newTheta[3]);
//			//newTheta[2]=max(0.0f, (Nmax-newTheta[3])*2*pi*PSFSigma_x*PSFSigma_y*sqrt(2.0f));
//			//newTheta[4]=0;
//
//		case 4:
//			//float newClamp[NV]={1.0,1.0,100,20,0.5,0.5}, oldClamp[NV]={1.0,1.0,100,20,0.5,0.5};
//
//			//kernel_CenterofMass2D(sz, s_data, &newTheta[0], &newTheta[1]);
//			//kernel_GaussFMaxMin2D(sz, PSFSigma, s_data, &Nmax, &newTheta[3]);
//			//newTheta[2]=max(0.0f, (Nmax-newTheta[3])*2*pi*PSFSigma*PSFSigma);
//			//newTheta[4]=PSFSigma;
//			//newTheta[5]=PSFSigma;
//
//		case 5:
//			//float newClamp[NV]={1.0,1.0,100,20,2.0}, oldClamp[NV]={1.0,1.0,100,20,2.0};
//
//			//kernel_CenterofMass2D(sz, s_data, &newTheta[0], &newTheta[1]);
//			//kernel_GaussFMaxMin2D(sz, 1.5, s_data, &Nmax, &newTheta[3]);
//			//newTheta[2]=max(0.0, (Nmax-newTheta[3])*2*pi*1.5*1.5);
//			//newTheta[4]=float(spline_zsize)/2;
//
//		}
//
//		float M[NV*NV],Diag[NV], Minv[NV*NV];
//		float newSign[NV] = {0}, oldSign[NV] = {0};
//		float newUpdate[NV] = {0},oldUpdate[NV] = {0};	
//		float newDudt[NV] ={0};
//		//float *M=(float*)malloc(NV*NV*sizeof(float));
//		//float *Diag=(float*)malloc(NV*sizeof(float));
//		//float *Minv=(float*)malloc(NV*NV*sizeof(float));
//
//		//float *newSign=(float*)malloc(NV*sizeof(float));
//		//float *oldSign=(float*)malloc(NV*sizeof(float));
//		//float *newUpdate=(float*)malloc(NV*sizeof(float));
//		//float *oldUpdate=(float*)malloc(NV*sizeof(float));
//		//float *newDudt=(float*)malloc(NV*sizeof(float));
//		//float *jacobian=(float*)malloc(NV*sizeof(float));
//		//float *hessian=(float*)malloc(NV*NV*sizeof(float));
//		//float *L=(float*)malloc(NV*NV*sizeof(float));
//		//float *U=(float*)malloc(NV*NV*sizeof(float));
//
//		//memset(M,0,NV*NV*sizeof(float));
//		//memset(Diag,0,NV*sizeof(float));
//		//memset(Minv,0,NV*NV*sizeof(float));
//		//memset(newSign,0,NV*sizeof(float));
//		//memset(oldSign,0,NV*sizeof(float));
//		//memset(newUpdate,0,NV*sizeof(float));
//		//memset(oldUpdate,0,NV*sizeof(float));
//		//memset(newDudt,0,NV*sizeof(float));
//		//memset(jacobian,0,NV*sizeof(float));
//		//memset(hessian,0,NV*NV*sizeof(float));
//		//memset(L,0,NV*NV*sizeof(float));
//		//memset(U,0,NV*NV*sizeof(float));
//
//		float model, cf, df, data;
//		float Div;
//		float PSFy, PSFx;
//		
//		float newLambda = 1.0, oldLambda = 1.0;
//
//		float newErr = 1e12, oldErr = 1e13;
//
//		
//		float jacobian[NV]={0};
//		float hessian[NV*NV]={0};
//		float t1,t2;
//		
//		float temp;
//		
//		int info;
//		//float L[NV*NV] = {0}, U[NV*NV] = {0};
//			float *L= new float [NV*NV];
//		float *U=new float [NV*NV];
//				cudaMemset (L,0,NV*NV*sizeof(float));
//		cudaMemset (U,0,NV*NV*sizeof(float));
//
//		
//
//		for (ii=0;ii<NV*NV;ii++)M[ii]=0;
//		for (ii=0;ii<NV*NV;ii++)Minv[ii]=0;
//		for (ii=0;ii<NV;ii++)oldTheta[ii]=newTheta[ii];
//		
//
//		//updateFitValues3D
//		newErr = 0;
//		memset(jacobian,0,NV*sizeof(float));
//		memset(hessian,0,NV*NV*sizeof(float));
//		for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {            
//			
//			kernel_DerivativeGauss2D_sigma(ii,jj,newTheta,newDudt,&model);
//			data=s_data[sz*jj+ii];
//
//			if (data>0)
//				newErr = newErr + 2*((model-data)-data*log(model/data));
//			else if (data==0)
//				newErr = newErr + 2*model;
//
//			t1 = 1-data/model;
//			for (l=0;l<NV;l++){
//				jacobian[l]+=t1*newDudt[l];
//			}
//
//			t2 = data/pow(model,2);
//			for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
//				hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
//				hessian[m*NV+l] = hessian[l*NV+m];
//			}
//		}
//		//addPeak
//
//		//copyFitData
//
//		for (kk=0;kk<iterations;kk++) {//main iterative loop
//
//			if(abs((newErr-oldErr)/newErr)<TOLERANCE){
//				//newStatus = CONVERGED;
//				break;
//			}
//			else{
//				if(newErr>1.5*oldErr){
//					//copy Fitdata
//
//					for (i=0;i<NV;i++){
//						newSign[i]=oldSign[i];
//						newClamp[i]=oldClamp[i];
//						newTheta[i]=oldTheta[i];
//					}
//					newLambda = oldLambda;
//					newErr = oldErr;
//
//					newLambda = 10*newLambda;
//				}
//				else if(newErr<oldErr){
//					if (newLambda>1){
//						newLambda = newLambda*0.8;
//					}
//					else if(newLambda<1){
//						newLambda = 1;
//					}
//				}
//
//
//				for (i=0;i<NV;i++){
//					hessian[i*NV+i]=hessian[i*NV+i]*newLambda;
//				}
//				memset(L,0,NV*sizeof(float));
//				memset(U,0,NV*sizeof(float));
//				info = kernel_cholesky(hessian,NV,L,U);
//				if (info ==0){
//					kernel_luEvaluate(L,U,jacobian,NV,newUpdate);
//					//copyFitData
//					for (i=0;i<NV;i++){
//						oldSign[i]=newSign[i];
//						oldClamp[i]=newClamp[i];
//
//						oldTheta[i]=newTheta[i];
//					}
//					oldLambda = newLambda;
//					oldErr=newErr;
//
//
//					//updatePeakParameters
//					for (ll=0;ll<NV;ll++){
//						if (newSign[ll]!=0){
//							if (newSign[ll]==1&&newUpdate[ll]<0){
//								newClamp[ll]=newClamp[ll]*0.5;
//							}
//							else if (newSign[ll]==-1&&newUpdate[ll]>0){
//								newClamp[ll] = newClamp[ll]*0.5;
//							}
//						}
//
//						if (newUpdate[ll]>0){
//							newSign[ll]=1;
//						}
//						else{
//							newSign[ll]=-1;
//						}
//
//						newTheta[ll] = newTheta[ll]-newUpdate[ll]/(1+abs(newUpdate[ll]/newClamp[ll]));
//					}
//
//					newTheta[0] = max(newTheta[0],(float(sz)-1)/2-2);
//					newTheta[0] = min(newTheta[0],(float(sz)-1)/2+2);
//					newTheta[1] = max(newTheta[1],(float(sz)-1)/2-2);
//					newTheta[1] = min(newTheta[1],(float(sz)-1)/2+2);
//					newTheta[2] = max(newTheta[2],1.0);
//					newTheta[3] = max(newTheta[3],0.01);
//					newTheta[4] = max(newTheta[4],0.0);
//					newTheta[4] = min(newTheta[4],sz/2.0f);
//
//
//					newErr = 0;
//					memset(jacobian,0,NV*sizeof(float));
//					memset(hessian,0,NV*NV*sizeof(float));
//					for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
//						//PSFx=kernel_IntGauss1D(ii, newTheta[0], newTheta[4]);
//						//PSFy=kernel_IntGauss1D(jj, newTheta[1], newTheta[4]);
//
//						//model=newTheta[3]+newTheta[2]*PSFx*PSFy;
//						//data=s_data[sz*jj+ii];
//
//						////calculating derivatives
//						//kernel_DerivativeIntGauss1D(ii, newTheta[0], newTheta[4], newTheta[2], PSFy, &newDudt[0],NULL);
//						//kernel_DerivativeIntGauss1D(jj, newTheta[1], newTheta[4], newTheta[2], PSFx, &newDudt[1],NULL);
//						//kernel_DerivativeIntGauss2DSigma(ii, jj, newTheta[0], newTheta[1], newTheta[4], newTheta[2], PSFx, PSFy, &newDudt[4],NULL);
//						//newDudt[2] = PSFx*PSFy;
//						//newDudt[3] = 1.0f;
//						/*(*jacf)(ii,jj,newTheta,newDudt,&model);*/
//						kernel_DerivativeGauss2D_sigma(ii,jj,newTheta,newDudt,&model);
//						//(*jacf)(ii,jj,newTheta,newDudt,&model);
//						data=s_data[sz*jj+ii];
//
//						if (data>0)
//							newErr = newErr + 2*((model-data)-data*log(model/data));
//						else if (data==0)
//							newErr = newErr + 2*model;
//
//						t1 = 1-data/model;
//						for (l=0;l<NV;l++){
//							jacobian[l]+=t1*newDudt[l];
//						}
//
//						t2 = data/pow(model,2);
//						for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
//							hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
//							hessian[m*NV+l] = hessian[l*NV+m];
//						}
//					}
//				}
//				else
//				{
//					newLambda = 10*newLambda;
//				}
//
//				//copyFitdata
//			}
//		}
//
//		// Calculating the CRLB and LogLikelihood
//		Div=0.0f;
//		for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
//			//need to check why don't use newTheta[4] instead of PSFSigma!!!
//			PSFx=kernel_IntGauss1D(ii, newTheta[0], PSFSigma);
//			PSFy=kernel_IntGauss1D(jj, newTheta[1], PSFSigma);
//
//			model=newTheta[3]+newTheta[2]*PSFx*PSFy;
//			data=s_data[sz*jj+ii];
//
//			//calculating derivatives
//			kernel_DerivativeIntGauss1D(ii, newTheta[0], newTheta[4], newTheta[2], PSFy, &newDudt[0], NULL);
//			kernel_DerivativeIntGauss1D(jj, newTheta[1], newTheta[4], newTheta[2], PSFx, &newDudt[1], NULL);
//			kernel_DerivativeIntGauss2DSigma(ii, jj, newTheta[0], newTheta[1], newTheta[4], newTheta[2], PSFx, PSFy, &newDudt[4], NULL);
//			newDudt[2] = PSFx*PSFy;
//			newDudt[3] = 1.0f;
//
//			//Building the Fisher Information Matrix
//			for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
//				M[kk*NV+ll]+= newDudt[ll]*newDudt[kk]/model;
//				M[ll*NV+kk]=M[kk*NV+ll];
//			}
//
//			//LogLikelyhood
//			if (model>0)
//				if (data>0)Div+=data*log(model)-model-data*log(data)+data;
//				else
//					Div+=-model;
//		}
//
//		// Matrix inverse (CRLB=F^-1) and output assigments
//		kernel_MatInvN(M, Minv, Diag, NV);
//		//write to global arrays
//		for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+BlockSize*bx+tx]=newTheta[kk];
//		for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+BlockSize*bx+tx]=Diag[kk];
//		d_LogLikelihood[BlockSize*bx+tx] = Div;
//
//		return;
//}
//
