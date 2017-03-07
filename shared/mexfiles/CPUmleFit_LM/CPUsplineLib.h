#ifndef CPUSPLINELIB_H
#define CPUSPLINELIB_H

//__device__ void kernel_computeDelta3D(float z_delta, float y_delta, float x_delta, float *delta_f, float *delta_dxf, float *delta_ddxf, float *delta_dyf, float *delta_ddyf, float *delta_dzf, float *delta_ddzf);


 void kernel_computeDelta3D(float z_delta, float y_delta, float x_delta, float *delta_f, float *delta_dxf, float *delta_dyf, float *delta_dzf);


 float kernal_fAt3D(int zc, int yc, int xc, int xsize, int ysize, int zsize, float *delta_f, const float *coeff);

 int kernel_cholesky(float *A,int n, float *L, float*U);

 void kernel_luEvaluate(float *L,float *U, float *b, int n, float *x);


 void kernel_DerivativeSpline(int zc, int yc, int xc, int xsize, int ysize, int zsize, float *delta_f, float *delta_dxf, float *delta_dyf, float *delta_dzf,const float *coeff,float *theta, float*dudt,float*model);


#endif