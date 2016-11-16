#ifndef GPUSPLINELIB_H
#define GPUSPLINELIB_H

//__device__ void kernel_computeDelta3D(float z_delta, float y_delta, float x_delta, float *delta_f, float *delta_dxf, float *delta_ddxf, float *delta_dyf, float *delta_ddyf, float *delta_dzf, float *delta_ddzf);


__device__ void kernel_computeDelta3D(float z_delta, float y_delta, float x_delta, float *delta_f, float *delta_dxf, float *delta_dyf, float *delta_dzf);


__device__ float kernal_fAt3D(int zc, int yc, int xc, int xsize, int ysize, int zsize, float *delta_f, float *coeff);

__device__ int kernel_cholesky(float *A,int n, float *L, float*U);

__device__ void kernel_luEvaluate(float *L,float *U, float *b, int n, float *x);

#endif