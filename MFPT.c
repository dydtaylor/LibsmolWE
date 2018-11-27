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

#define NBINSMAX 1000
#define BINCONTENTSMAXMAX 2000
#define ISIMMAXMAX 1000000 

struct paramsOUWeightedEnsemble{
	unsigned int tau; //In integer units of dt: e.g. if dt=.01, tau=.1, then this value should be .1/.01 = 10;
	unsigned int repsPerBin; // Target number of replicas per bin. Sometimes called "MTarg".
	unsigned int tauMax; //Number of weighted ensemble steps to do, max sim time.
			// If tauMax = 1000, tau = 10, dt = .005, then sim will end at t=tauMax*tau*dt = 50;
	unsigned int nBins;
	unsigned int fluxBin;
	double binDefs[NBINSMAX];
};

struct paramsOUDynamicsEngine{
	// Simulation parameters
	double dt; // Timestep of dynamics engine
	double worldLength;
	double roiR;
	double difC;
	int nPart;
	// Model parameters
	double tauSlow; // Timescale of OU process
	double sigmaX; // Standard deviation at steady-state of OU process
};

struct replicas{
	simptr sims[ISIMMAXMAX];
	double weights[ISIMMAXMAX];
	unsigned int binLocs[ISIMMAXMAX];
	unsigned int nBins;
	unsigned int binContentsMax[NBINSMAX];
	unsigned int binContents[BINCONTENTSMAXMAX][NBINSMAX];
	unsigned int iSimMax;
};

struct paramsOUWeightedEnsemble paramsWeOu;
struct paramsOUDynamicsEngine paramsDeOu;
struct replicas Reps;
FILE *errFile;

void initialDist(int nInit){
	int jSim;
	double lowBounds[] = {-paramsDeOu.worldLength/2,-paramsDeOu.worldLength/2};
	double highBounds[] = {paramsDeOu.worldLength/2, paramsDeOu.worldLength/2};
	double botLeftCornerRect[] = {-paramsDeOu.worldLength/2, -paramsDeOu.worldLength/2, paramsDeOu.worldLength};
	double topRightCornerRect[] = {paramsDeOu.worldLength/2, paramsDeOu.worldLength/2, -paramsDeOu.worldLength};
	double roiParams[] = {0.0, 0.0, 1, 30};
	double insideRoi[] = {0.0, 0.0};
	for(jSim = 0; jSim < nInit; jSim++){
		//Molecules + BCs
		Reps.sims[jSim] = smolNewSim(2, lowBounds,highBounds);
		smolSetGraphicsParams(Reps.sims[jSim], "none", 1, 0);
		smolSetSimTimes(Reps.sims[jSim],0,10000,paramsDeOu.dt);
		smolAddSpecies(Reps.sims[jSim],"A",NULL);
		smolSetSpeciesMobility(Reps.sims[jSim],"A",MSall, paramsDeOu.difC, 0, 0);
		smolAddSolutionMolecules(Reps.sims[jSim], "A", paramsDeOu.nPart, NULL, NULL);
		smolAddPanel(Reps.sims[jSim], "bounds", PSrect, NULL, "-x", topRightCornerRect);
		smolAddPanel(Reps.sims[jSim], "bounds", PSrect, NULL, "-y", botLeftCornerRect);
		smolAddPanel(Reps.sims[jSim], "bounds", PSrect, NULL, "+x", botLeftCornerRect);
		smolAddPanel(Reps.sims[jSim], "bounds", PSrect, NULL, "+y", topRightCornerRect);
		smolSetSurfaceAction(Reps.sims[jSim], "bounds", PFboth, "all", MSall, SAreflect);
		
		//ROI Surface + compartment
		smolAddSurface(Reps.sims[jSim], "roi");
		smolAddPanel(Reps.sims[jSim],"roi", PSsph, NULL, "+0", roiParams);
		smolSetSurfaceAction(Reps.sims[jSim], "roi", PFboth, "all", MSall, SAtrans);
		smolAddCompartment(Reps.sims[jSim],"roiComp");
		smolAddCompartmentSurface(Reps.sims[jSim],"roiComp","roi");
		smolAddCompartmentPoint(Reps.sims[jSim],"roiComp",insideRoi);
		
		smolAddCommandFromString(Reps.sims[jSim], "e ifincmpt A = 0 roiComp stop");
		smolUpdateSim(Reps.sims[jSim]);
	}
}

void dynamicsEngine(){
	int jSim;
	double currentT, breakT;
	for(jSim = 0; jSim < Reps.iSimMax; jSim++){
		currentT = Reps.sims[jSim]->time;
		breakT = currentT + paramsWeOu.tau*paramsDeOu.dt;
		smolRunSimUntil(Reps.sims[jSim], breakT);
	}
}

int findBin(){
	int nInBin;
	double molX, molY;
	nInBin = 0;
	for(nMol = 0; nMol < paramsDeOu.nPart; nMol++){
		molX = Reps.sims[jSim]->mols->live[0][nMol]->pos[0]; //first index changes based on how we org mols
		molY = Reps.sims[jSim]->mols->live[0][nMol]->pos[1];
		if(molX^2 + molY^2 < paramsDeOu.roiR){ //Needs to be changed to deal with varying ROIs
			nInBin++;
		}
	}
	return nInBin;
}

void splitMerge(){
	
}

double fluxes(){
	
}

int main(int argc, char *argv[]){
	int tauQuarter, tauMax, rngBit, iBin, nWE;
	
	FILE *DEFile, *WEFile, *BINFile, *FLFile, *SIMFile;
	DEFile = fopen("dynamicsParams.txt","r");
	WEFile = fopen("WEParams.txt","r");
	BINFile = fopen("Bins.txt","r");
	errFile = fopen(argv[3], "w");
	fclose(DEFile);
	fclose(WEFile);
	fclose(BINFile);
	fclose(errFile);
	return 0;
}