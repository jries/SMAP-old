#include <windows.h>
#pragma comment(lib, "kernel32.lib")

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>
#include <thread>
#include <mutex>
using namespace std;

#ifndef max
//! not defined in the C standard used by visual studio
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
//! not defined in the C standard used by visual studio
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif
#define pi 3.141592f


std::mutex worker_mutex;

void worker(const size_t id, const size_t NFitraw, size_t &current) {
	
	size_t subregion=0; 

	worker_mutex.lock();
	while (current < NFitraw) {
		subregion = current++;
		worker_mutex.unlock();

		mexPrintf("thread %d processing subregion %d\n",id,subregion);
		Sleep(50);
		worker_mutex.lock();
	}
	worker_mutex.unlock();
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
	//declare all vars

	//check for required inputs, correct types, and dimensions
	//1D vectors still return 2D

	//retrieve all inputs

	//validate input values(this section better not be blank!)

	//do stuff
	const size_t Nfitsraw = 100;
	size_t current = 0;

	SYSTEM_INFO sysinfo;
	GetSystemInfo( &sysinfo );

	const int numCPU = sysinfo.dwNumberOfProcessors;
	const int threadlimit = std::thread::hardware_concurrency();

	std::thread *t = new std::thread[numCPU];
	for (size_t ii=0; ii<numCPU; ii++)
		t[ii] = std::thread(worker, ii, Nfitsraw, std::ref(current));
	
	for (size_t ii=0; ii<numCPU; ii++)
		t[ii].join();

	mexPrintf("cores %d thread limit %d\n", numCPU,threadlimit);
	mexPrintf("%d regions out of %d regions processed\n",current, Nfitsraw);

	//allocate results memory and copy results back to matlab


	//clean up anything not using an mxMalloc

	return;
 }