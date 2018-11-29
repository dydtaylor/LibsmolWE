#include "weSmoldyn.h"
#include "smolDynamics.c"

void splitMerge(){
	int binMin = 0;
	int binMax = Reps.nBins;
	unsigned int dummyInd;/*Just a dummy index that will hold the current index location in both split and merge loop*/
	int rowCol[3]; /* YET ANOTHER dummy that tells me which columns in the bincontents row will be combined together, with the third element giving the deleted column*/
	int mergeInd[2]; /* Array containing indices of 2 elements to be merged*/
	int keptInd[2]; /*Previous array, reordered s.t. the first index is the one that will be kept*/
	int splitInd; /* Index of replica to split*/
	int splitBin, mergeBin, repInBin, entryCheck1, entryCheck2, iSimMaxReplace; //Looping variables

	double p0; /*Doubles giving the merge probabilities*/
	double randPull; /*Double storing a pull from RAND*/
	//unsigned int binConPrev[BINCONTENTSMAXMAX][NBINSMAX];

	/*Merging loop*/
	for(mergeBin = binMin; mergeBin <binMax; mergeBin++){
		/*Initialize the merge indices with the first 2 indices stored in appropriate BinContents row*/

		while((Reps.binContentsMax[mergeBin]>paramsWe.repsPerBin) && Reps.binContentsMax[mergeBin] > 0){
			mergeInd[0] = Reps.binContents[0][mergeBin];
			mergeInd[1] = Reps.binContents[1][mergeBin];
			rowCol[0] = 0;
			rowCol[1] = 1;
			/*Find the locations of the two smallest weights and combine them together*/
			for(repInBin = 0; repInBin < Reps.binContentsMax[mergeBin];repInBin++){ //NUM_Rows needs to change to an appropriate BCM statement

				dummyInd = Reps.binContents[repInBin][mergeBin];
				/*If the weight of this index is greater than the weight in the first merge index,
				but smaller than the weight in the second index, then replace the first merge index with the new index

				Otherwise, if the weight of the dummy index is greater than the weight in the 2nd index, then replace the second
				index with the new index
				*/
				if((Reps.weights[dummyInd] < Reps.weights[mergeInd[0]])&&(Reps.weights[dummyInd] > Reps.weights[mergeInd[1]])){
					mergeInd[0] = dummyInd;
					rowCol[0] = repInBin;

				}
				else if((Reps.weights[dummyInd] < Reps.weights[mergeInd[1]])&& mergeInd[0] != dummyInd){
					mergeInd[1] = dummyInd;
					rowCol[1] = repInBin;
				}
			}

			if(mergeInd[0] == mergeInd[1]){
					//fprintf(errFile, "Merge Error \n");
			}

			/*Decide which index to keep*/
			p0 = Reps.weights[mergeInd[0]] / (Reps.weights[mergeInd[0]]+Reps.weights[mergeInd[1]]);
			randPull = RAND;
			if(randPull<p0){ // JUN COMMENT: Do you need a variable randPull? Just use RAND, and, below, use else instead of else-if.
				keptInd[0] = mergeInd[0];
				keptInd[1] = mergeInd[1];
				rowCol[2] = rowCol[1];
			}
			else if(randPull > p0){
				keptInd[0] = mergeInd[1];
				keptInd[1] = mergeInd[0];
				rowCol[2] = rowCol[0];
			}

			/*Update weight of the kept index*/
			if( (isnan(Reps.weights[keptInd[0]])) || (isnan(Reps.weights[keptInd[1]]))){ // JUN COMMENT: Is it possible to use isnan()?
				//fprintf(errFile, "WARNING: Moving NAN Weight \n");
			}
			Reps.weights[keptInd[0]] = Reps.weights[keptInd[0]] + Reps.weights[keptInd[1]];

			/*Replace the old simulation with the final non-NAN simulation*/

			if(isnan(Reps.weights[Reps.iSimMax])){
				//fprintf(errFile, "WARNING: Moving NAN Weight \n");
			}

			int simMaxHolder;

			/*Find iSimMax in binContents and replace it with keptInd[1]*/
			/*For loop: iterates through the row of the bin ISM is in.
			If: Checks to see if the index of the row we're iterating through is ISM. If it is, it replaces it with KI[1], the deleted index.

			The reason we replace ISM with the deleted index is for the incrementing nature of keeping track of the sims.
			We aren't deleting ISM, we're moving it to the spot where the deleted sim was.

			This has problems when the deleted sim is ISM. When that happens, I will take the last element of the table of contents and
			move it to where ISM is in the table.*/

			for(iSimMaxReplace = 0; iSimMaxReplace < Reps.binContentsMax[Reps.binLocs[Reps.iSimMax]];iSimMaxReplace++){
				if(Reps.binContents[iSimMaxReplace][Reps.binLocs[Reps.iSimMax]] == Reps.iSimMax){
					Reps.binContents[iSimMaxReplace][Reps.binLocs[Reps.iSimMax]] = keptInd[1];
					simMaxHolder = iSimMaxReplace;
					iSimMaxReplace = Reps.binContentsMax[Reps.binLocs[Reps.iSimMax]];
				}
			}

			if (Reps.iSimMax ==keptInd[1]){
				Reps.binContents[simMaxHolder][Reps.binLocs[Reps.iSimMax]] = Reps.binContentsMax[Reps.binLocs[Reps.iSimMax]] - 1;
			}

			Reps.sims[keptInd[1]] = Reps.sims[Reps.iSimMax];
			Reps.weights[keptInd[1]] = Reps.weights[Reps.iSimMax];
			Reps.binLocs[keptInd[1]] = Reps.binLocs[Reps.iSimMax];

			//Setting things to NAN is technically unnecessary
			/*Remove the duplicate non-NAN simulation at the end of the non-NANs
			Reps.sims[Reps.iSimMax] = NAN;*/
			Reps.weights[Reps.iSimMax] = NAN;
			Reps.binLocs[Reps.iSimMax] = NAN;
			Reps.iSimMax--;

			/*Reorganize the binContents matrix*/
			Reps.binContents[rowCol[2]][mergeBin] = Reps.binContents[Reps.binContentsMax[mergeBin] -1][mergeBin];
			Reps.binContents[-1+Reps.binContentsMax[mergeBin]][mergeBin] = NAN;
			Reps.binContentsMax[mergeBin]--;

			for(entryCheck1 = 0; entryCheck1 < Reps.binContentsMax[mergeBin]-1;entryCheck1++){
				for(entryCheck2 = entryCheck1 + 1; entryCheck2< Reps.binContentsMax[mergeBin]; entryCheck2++){
					if(Reps.binContents[entryCheck1][mergeBin]==Reps.binContents[entryCheck2][mergeBin]){
						//fprintf(errFile, "ERROR: Duplicate entries in BC \n");
					}
				} // finished j loop through this bin's contents
			} // finished i loop through this bin's contents
		} // finished while-loop to check this bin
	} // finished merging loop through bins

	/*Splitting Loop*/
	for(splitBin = binMin; splitBin < binMax; splitBin++){
		while((Reps.binContentsMax[splitBin]<paramsWe.repsPerBin)&&(Reps.binContentsMax[splitBin]>0)){
			splitInd = Reps.binContents[0][splitBin];
			for(repInBin = 0; repInBin < Reps.binContentsMax[splitBin];repInBin++){
				dummyInd = Reps.binContents[repInBin][splitBin];
				if(Reps.weights[dummyInd]>Reps.weights[splitInd]){
					splitInd = dummyInd;
				}
			}
			Reps.sims[Reps.iSimMax+1] = Reps.sims[splitInd]; //Needs to create a new Sim, this just copies a pointer which will never work
			Reps.weights[splitInd] = Reps.weights[splitInd] / 2;
			Reps.weights[Reps.iSimMax+1] = Reps.weights[splitInd];
			Reps.binLocs[Reps.iSimMax+1] = Reps.binLocs[splitInd];
			Reps.iSimMax++;
			if(Reps.iSimMax > ISIMMAXMAX){
				printf("ERROR: iSimMax out of bounds");
			}
			Reps.binContents[-1+Reps.binContentsMax[splitBin]][splitBin] = Reps.iSimMax;
			Reps.binContentsMax[splitBin]++;
		}
	} // finished splitting loop
	return;
}


double fluxes(){
	double fluxOut = 0;
	double weightSum;
	int jWeight, iReps, iSimMaxReplace, iSim;
	int nFlux = Reps.binContentsMax[paramsWe.fluxBin]; // Number of replicas in the fluxBin

	weightSum = 0;
	for(jWeight = 0; jWeight <= Reps.iSimMax; jWeight++){
		weightSum = weightSum + Reps.weights[jWeight];
	}

	if(fabs(weightSum-1)>1e-6){
		printf("Weight not conserved \n");
	}
	
	// loop through replicas in flux bin and delete them
	if(nFlux>0){
		for(iReps = nFlux -1; iReps >=0; iReps--){
			 // JUN COMMENT: Is this a validation check? Shouldn't all these be in the fluxbin?
				fluxOut += Reps.weights[Reps.binContents[iReps][paramsWe.fluxBin]]; // JUN COMMENT: Use += notation?
				Reps.sims[Reps.binContents[iReps][paramsWe.fluxBin]] = Reps.sims[Reps.iSimMax];
				Reps.weights[Reps.binContents[iReps][paramsWe.fluxBin]] = Reps.weights[Reps.iSimMax];
				Reps.binLocs[Reps.binContents[iReps][paramsWe.fluxBin]] = Reps.binLocs[Reps.iSimMax];

				
				for(iSimMaxReplace = 0; iSimMaxReplace < Reps.binContentsMax[Reps.binLocs[Reps.iSimMax]];iSimMaxReplace++){
					if(Reps.binContents[iSimMaxReplace][Reps.binLocs[Reps.iSimMax]] == Reps.iSimMax){
						Reps.binContents[iSimMaxReplace][Reps.binLocs[Reps.iSimMax]] = Reps.binContents[iReps][paramsWe.fluxBin];
					}
				}

				//Reps.sims[Reps.iSimMax] = NAN;
				Reps.weights[Reps.iSimMax] = NAN;
				Reps.binLocs[Reps.iSimMax] = NAN;
				Reps.iSimMax--;
				Reps.binContentsMax[paramsWe.fluxBin]--;
		}
	}

	if(Reps.binContentsMax[paramsWe.fluxBin] != 0){
			printf("Non Zero weight in flux bin \n");
	};


	if(fluxOut != 0){

			for(iSim = 0; iSim < Reps.iSimMax; iSim++){
				Reps.weights[iSim] = Reps.weights[iSim]/(weightSum-fluxOut);
			}

	}

	return fluxOut;
}

int main(int argc, char *argv[]){
	int tauQuarter, tauMax, rngBit, iBin, nWE;
	
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

	printf("Parameters loaded\n");

	tauQuarter = paramsWe.tauMax / 4;
	tauMax = paramsWe.tauMax;

	printf("Tau loops + flux vector made \n");

	rngBit = atoi(argv[4]);
	fprintf(errFile, "iseed=%lx\n", RanInitReturnIseed(rngBit));
    fclose(errFile);
	initialDist(paramsWe.repsPerBin);

	printf("Initial Distribution Made \n");
	printf("Initial Bin Location = %i \n", Reps.binLocs[0]);
	
	for(nWE = 0; nWE < tauQuarter; nWE++){
		printf("Tau Step: %i \n", nWE);
		splitMerge();
		for(int iBin = 0; iBin < Reps.nBins; iBin++){
			Reps.binContentsMax[iBin] = 0;
		}
		for(int iSim = 0; iSim < Reps.iSimMax; iSim++){
			printf("\n \n \n sim = %i \n \n \n",iSim);
			dynamicsEngine(Reps.sims[iSim]);
			Reps.binLocs[iSim] = findBin(Reps.sims[iSim]);
			Reps.binContents[Reps.binContentsMax[Reps.binLocs[iSim]]][Reps.binLocs[iSim]] = iSim;
			Reps.binContentsMax[Reps.binLocs[iSim]]++;
		}
	}
	
	for(nWE = tauQuarter; nWE < tauMax; nWE++){
		printf("Tau Step: %i \n", nWE);
		FLFile = fopen(argv[2],"a");
		fprintf(FLFile, "%E \n", fluxes());
		fclose(FLFile);
		splitMerge();
		for(int iBin = 0; iBin < Reps.nBins; iBin++){
			Reps.binContentsMax[iBin] = 0;
		}
		for(int iSim = 0; iSim < Reps.iSimMax; iSim++){
			dynamicsEngine(Reps.sims[iSim]);
			Reps.binLocs[iSim] = findBin(Reps.sims[iSim]);
			Reps.binContents[Reps.binContentsMax[Reps.binLocs[iSim]]][Reps.binLocs[iSim]] = iSim;
			Reps.binContentsMax[Reps.binLocs[iSim]]++;
		}
	}
	
	
	SIMFile = fopen(argv[1], "w");
	for(iBin = 0;iBin <= Reps.iSimMax; iBin++){ // JUN COMMENT: Better name for counter, like iBin.
		fprintf(SIMFile, "%i, %E \n", findBin(Reps.sims[iBin]), Reps.weights[iBin]);
	}
	fclose(SIMFile);
	
	
	return 0;
}