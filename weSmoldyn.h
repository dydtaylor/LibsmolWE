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

#define NBINSMAX 300
#define BINCONTENTSMAXMAX 50000
#define ISIMMAXMAX 50000
#define NULLDEVICE "/dev/null"

#define SMOLTIMEMAX 10000

#define NFLUXBINS 150
#define KSCRITICAL .02
#define DKSCRITICAL 0.3
#define DEBUGGING 0 //0 disables debugging text output.
#define STOPCOMMAND 1
#define WEENABLE 1 //Ignores WE and only runs Smoldyn dynamics
#define ROBINS 1 // Specifies WE binning is based on N inside the ROI order parameter
#define MONOFRACEACHDT 0 //Determines frequency of measuring monomerization fraction
#define FIXEDTIME 1
#define NNZMIN 500 //DKS Test non-zero min for 

struct paramsWeightedEnsemble{
	unsigned int tau; //In integer units of dt: e.g. if dt=.01, tau=.1, then this value should be .1/.01 = 10;
	unsigned int repsPerBin; // Target number of replicas per bin. Sometimes called "MTarg".
	unsigned int tauMax; //Number of weighted ensemble steps to do, max sim time.
			// If tauMax = 1000, tau = 10, dt = .005, then sim will end at t=tauMax*tau*dt = 50;
	unsigned int nBins;
	int fluxBin;
	unsigned int nInit;
	//double binDefs[NBINSMAX];
} ;

struct paramsDynamicsEngine{
	// Simulation parameters
	double dt; // Timestep of dynamics engine
	double worldLength;
	double roiR;
	double difM;
	double difD;
	double bindR;
	double unbindK;
	int nPart;
	int reactBit;
	int reentryRateBit;
	double reentryRate;
	int corralsBit;
	double corralWidth;
	double corralRate;
	int densityBit;
	double density;
	int monomerStart;
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

struct fluxCounts{
	double fluxMax;
	double binDefs[NFLUXBINS+1];
	double ksStat;
	double dualKS;
	int binCounts[NFLUXBINS];
	int oldCounts[NFLUXBINS];
	int dualCounts[NFLUXBINS];
	int oldDualCounts[NFLUXBINS];
	int nT;
};

struct BinDefinitions{
	double binDefArray[4*NBINSMAX];
	int nBins;
	int currentDims; //Only implemented to be 1 right now
	int customBins;
};

struct paramsWeightedEnsemble paramsWe;
struct paramsDynamicsEngine paramsDe;
struct replicas Reps;
struct fluxCounts fluxCDF;
struct BinDefinitions binDefs;

double lowBounds[2];
double highBounds[2];
double botLeftCornerRect[3];
double topRightCornerRect[3];
double roiParams[4];
double insideRoi[2];

/*Document Explanation (At bottom to prevent interfering with sed commands)
!------------------------------------------------------------------------------
!This header file defines several structures as well as macros for executing
!the program, as well as declare several global variables. As all this does
!is define these objects, I will describe them each here.
!Macros:
!	RAND, RANDINT: Calls to twister.c random numbers, RAND pulls from U[0,1]
!		and RANDINT from 32bit integers
!	PI: Circumference/diameter for a circle
!	NBINSMAX, BINCONTENTSMAXMAX,ISIMMAXMAX: Macros defined to determine the 
!		size of the replicas struct. The replicas struct would ideally change
!		based on parameters of the system, but C doesn't allow that to be done
!		easily, so we set "large enough" values for these parameters. Detailed
!		explanation on where these sizes come from will be listed in the struct
!		definition
!	NULLDEVICE: For file output, in case any output wants to be suppressed.
!		Currently not used.
!	SMOLTIMEMAX: Smoldyn needs a maximum simulation time defined. This shouldn't
!		ever be reached in this program, so it's again set to be "large enough"
!		given the parameters I'm using for the system. Could conceivably be 
!		moved to a struct and loaded from a text file.
!	NFLUXBINS: Similar to NBINSMAX, arbitrary: large enough: value to be used
!		in the fluxCounts struct.
!	KSCRITICAL, DKSCRITICAL: Accuracy limits for the KS/DKS test in weSmoldyn.c
!		Also could be moved to a text file.
!	DEBUGGING: Enables more verbose output when set to 1.
!	STOPCOMMAND: Determines whether smoldyn simulations will terminate in between
!		weighted ensemble measurements. Should be 1, was included for debugging
!	WEENABLE: Enables/disables the weighted ensemble steps. Should be 1,
!		included for debugging.
!	ROBINS: When set to 1 uses the number of molecules inside the RoI as the
!		order parameter, when set to 0 uses the monomerization fraction.
!	MONOFRACEACHDT: Specifies how frequently monomerization fraction should be 
!		measured. Highly recommended to be 0 to avoid excessive output.
!	FIXEDTIME: Specifies whether or not to stop at tauMax (WEParams.txt) or
!		whether to stop once KS/DKS requirements have been met.
!	NNZMIN: Minimum number of non-zero measurements to pass the KS/DKS requirements
!Structs:
!	paramsWeightedEnsemble: Struct containing key parameters that are used for 
!		weighted ensemble steps. paramsWe is the only instance of this struct
!		tau: Period of weighted ensemble steps, in units of dynamics time steps
!		repsPerBin: How many replicas to be in occupied bins after WE reweighting
!		tauMax: How many weighted ensemble steps to execute.
!		nBins: How many bins to actually use for the order parameter
!		fluxBin: Which bin to measure the probability flux into
!		ninit: How many unique replicas to initialize for the system
!	paramsDynamicsEngine: Struct containing key dynamics engine parameters.
!		paramsDe is the only instance of this struct.
!		dt: How many natural time units for the system smoldyn should use for the 
!		its timesteps.
!		worldLength: How many natural distance units the domain will be long/wide
!		roiR: How many natural distance units the RoI radius will be.
!		difM: Monomer diffusion constant (in units of nat. dist^2 / nat time)
!		difD: Dimer diffusion constant (nat dist^2/nat time)
!		bindR: Binding radius for monomers (nat distance)
!		unbindK: Unbinding rate for dimers (1/nat time)
!		nPart: Maximum number of monomers in the system (twice as many max dimers)
!		reactBit: Enables binding / unbinding when set to 1, otherwise no rxns
!		reentryRateBit: Enables / disables variable reentry rates.
!		reentryRate: Rate for transmission through surface
!		densityBit, density: Enables/disables using density + nPart to set 
!		world length
!		monomerStart: Specifies whether simulation should start with monomers
!		or dimers
!	replicas: Contains information for indexing WE replicas, including weights,
!		bin locations, and total number of replicas currently
!		sims, weights, binLocs: Arrays of simpointers (pointers to an 
!		individual smoldyn simulation i.e. replica), replica weights, and
!		replica bin locations. weights[i] and binLocs[i] are each associated
!		with sims[i]
!		nBins: How many bins to use for order parameters
!		binContentsMax,binContents: "Table of contents" matrices and indices
!		binContents is a matrix tracks which replica is in which bin (e.g. sims[0] 
!		in 	bin 4, sims[3] in bin 5, etc). binContentsMax tells us how occupied
!		each column of the binContents matrix is. The sizes are set by macros
!		iSimMax: Tells us how many sims are currently being used / the last used 
!		index of the sims/weights/binLocs. This and BinContentsMax should be 
!		multiples of the repsPerBin WE parameter in between WE steps but will
!		change during Split/Merge and Flux measurement
!	fluxCounts: Struct for sorting flux measurements into bins for KS testing
!		fluxMax: Maximum measured flux (used for defining histogram bins)
!		binDefs: Definition of histogram bins used for KS testing
!		ksStat, dualKS: (daul) KS stat measured during the KS test
!		binCounts, oldCounts, dualCounts, oldDualCounts: Histogram counts
!		for flux measurements for each KS test. "old" counts has data from
!		2nd third measurements, other counts have data from 3rd third
!		nT: the tau step to perform the KS test at (integer)
!
!------------------------------------------------------------------------------
*/
