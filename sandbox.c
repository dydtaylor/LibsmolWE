#include "weSmoldyn.h"
#include "smolDynamics.c"
#include "smolSimCopy.c"

int main(int argc, char *argv[]){
	int tauQuarter, tauMax, rngBit, iBin, nWE, iSim;
	
	FILE *DEFile, *WEFile, *FLFile, *SIMFile, *errFile;
	DEFile = fopen("dynamicsParams.txt","r");
	WEFile = fopen("WEParams.txt","r");
	errFile = fopen(argv[3], "w");

	fscanf(WEFile,"%i %i %i %i", &paramsWe.tau, &paramsWe.repsPerBin, &paramsWe.tauMax, &paramsWe.nBins);
	fscanf(DEFile, "%lf %lf %lf %lf %i", &paramsDe.dt, &paramsDe.worldLength, &paramsDe.roiR, &paramsDe.difC, &paramsDe.nPart);
	paramsWe.fluxBin = 0;
	Reps.iSimMax = paramsWe.repsPerBin;
	fclose(DEFile);
	fclose(WEFile);
	Reps.nBins = paramsWe.nBins;

	printf("Parameters loaded\n");

	tauQuarter = paramsWe.tauMax / 4;
	tauMax = paramsWe.tauMax;

	printf("Tau loops + flux vector made \n");

	rngBit = atoi(argv[4]);
	fprintf(errFile, "iseed=%lx\n", RanInitReturnIseed(rngBit));
    fclose(errFile);
	initialDist(paramsWe.repsPerBin);
	for(iSim = 0; iSim < paramsWe.repsPerBin; iSim++){
		Reps.weights[iSim] = (double)1/paramsWe.repsPerBin;
		Reps.binLocs[iSim] = findBin(Reps.sims[iSim]);
		Reps.binContents[Reps.binContentsMax[Reps.binLocs[iSim]]][Reps.binLocs[iSim]] = iSim;
		Reps.binContentsMax[Reps.binLocs[iSim]]++;
	}
	
	for(iSim = 0; iSim < Reps.iSimMax; iSim++){
			printf("\n \n \n sim = %i \n \n \n",iSim);
			dynamicsEngine(Reps.sims[iSim]);
			Reps.binLocs[iSim] = findBin(Reps.sims[iSim]);
			Reps.binContents[Reps.binContentsMax[Reps.binLocs[iSim]]][Reps.binLocs[iSim]] = iSim;
			Reps.binContentsMax[Reps.binLocs[iSim]]++;
		}

	
	return 0
}