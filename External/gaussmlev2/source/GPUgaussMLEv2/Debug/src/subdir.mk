################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/mexFunction.cpp 

CU_SRCS += \
../src/GPUgaussMLEv2.cu \
../src/wrapper.cu 

CU_DEPS += \
./src/GPUgaussMLEv2.d \
./src/wrapper.d 

OBJS += \
./src/GPUgaussMLEv2.o \
./src/mexFunction.o \
./src/wrapper.o 

CPP_DEPS += \
./src/mexFunction.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	nvcc -I/Applications/MATLAB_R2010a.app/extern/include -G -g -O0 -m64 -gencode arch=compute_11,code=sm_11 -gencode arch=compute_20,code=sm_20 -gencode arch=compute_20,code=sm_21 -gencode arch=compute_30,code=sm_30 -odir "src" -M -o "$(@:%.o=%.d)" "$<"
	nvcc --compile -G -I/Applications/MATLAB_R2010a.app/extern/include -O0 -g -gencode arch=compute_11,code=compute_11 -gencode arch=compute_11,code=sm_11 -gencode arch=compute_20,code=compute_20 -gencode arch=compute_20,code=sm_20 -gencode arch=compute_20,code=sm_21 -gencode arch=compute_30,code=compute_30 -gencode arch=compute_30,code=sm_30 -m64  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	nvcc -I/Applications/MATLAB_R2010a.app/extern/include -G -g -O0 -m64 -gencode arch=compute_11,code=sm_11 -gencode arch=compute_20,code=sm_20 -gencode arch=compute_20,code=sm_21 -gencode arch=compute_30,code=sm_30 -odir "src" -M -o "$(@:%.o=%.d)" "$<"
	nvcc -I/Applications/MATLAB_R2010a.app/extern/include -G -g -O0 -m64 --compile  -x c++ -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


