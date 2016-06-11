#include <stdio.h>
#include <thread>
#include <mutex>


class FitInfo {
public:
	std::mutex *worker_mutex;
	volatile size_t current;
	size_t Nfitraw;
	size_t fittype;
	float* data;
	float PSFSigma;
	size_t sz;
	int iterations;
	float* Parameters;
	float* CRLBs;
	float* LogLikelihood;
	float Ax, Ay, Bx, By, gamma, d, PSFSigma_y;

	FitInfo() : 
		current(0),Nfitraw(0),data(0),PSFSigma(0),sz(0),iterations(10),Parameters(0),CRLBs(0),LogLikelihood(0) {
		worker_mutex = new std::mutex();
	};

	~FitInfo() {
		delete worker_mutex;
	}
};
