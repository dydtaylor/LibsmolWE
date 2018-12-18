#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <libsmoldyn.h>
#include <smoldyn.h>
#include "twister.c"

#define RAND genrand_real1()
#define RANDINT genrand_int32()
#define PI 3.14159265

#define NBINSMAX 100
#define BINCONTENTSMAXMAX 50000
#define ISIMMAXMAX 50000
#define NULLDEVICE "/dev/null"

struct paramsWeightedEnsemble{
	unsigned int tau; //In integer units of dt: e.g. if dt=.01, tau=.1, then this value should be .1/.01 = 10;
	unsigned int repsPerBin; // Target number of replicas per bin. Sometimes called "MTarg".
	unsigned int tauMax; //Number of weighted ensemble steps to do, max sim time.
			// If tauMax = 1000, tau = 10, dt = .005, then sim will end at t=tauMax*tau*dt = 50;
	unsigned int nBins;
	unsigned int fluxBin;
	unsigned int nInit;
	double binDefs[NBINSMAX];
} ;

struct paramsDynamicsEngine{
	// Simulation parameters
	double dt; // Timestep of dynamics engine
	double worldLength;
	double roiR;
	double difC;
	int nPart;
} ;

struct replicas{
	simptr sims[ISIMMAXMAX];
	double weights[ISIMMAXMAX];
	unsigned int binLocs[ISIMMAXMAX];
	unsigned int nBins;
	unsigned int binContentsMax[NBINSMAX];
	unsigned int binContents[BINCONTENTSMAXMAX][NBINSMAX];
	unsigned int iSimMax;
};

struct paramsWeightedEnsemble paramsWe;
struct paramsDynamicsEngine paramsDe;
struct replicas Reps;