#include "weSmoldyn.h"
#include "smolDynamics.c"

void splitMerge(int nWE){
	/*This function is largely unchanged from the OU Weighted Ensemble.
	Main difference is in the splitting loop, copySim was written to deal with smoldyn data structures
	*/
	int binMin = 0;
	int binMax = Reps.nBins;
	unsigned int dummyInd;/*Just a dummy index that will hold the current index location in both split and merge loop*/
	int rowCol[3]; /* YET ANOTHER dummy that tells me which columns in the bincontents row will be combined together, with the third element giving the deleted column*/
	int mergeInd[2]; /* Array containing indices of 2 elements to be merged*/
	int keptInd[2]; /*Previous array, reordered s.t. the first index is the one that will be kept*/
	int splitInd; /* Index of replica to split*/
	int splitBin, mergeBin, repInBin, entryCheck1, entryCheck2, iSimMaxReplace; //Looping variables

	int simMaxHolder;

	double p0; /*Doubles giving the merge probabilities*/
	double randPull; /*Double storing a pull from RAND*/
	//unsigned int binConPrev[BINCONTENTSMAXMAX][NBINSMAX];

	/*Merging loop*/
	for(mergeBin = binMin; mergeBin <binMax; mergeBin++){
		/*Initialize the merge indices with the first 2 indices stored in appropriate BinContents row*/

		while((Reps.binContentsMax[mergeBin]>paramsWe.repsPerBin) && Reps.binContentsMax[mergeBin] > 0){ // JUN COMMENT: What is the second condition for?
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

			if(mergeInd[0] == mergeInd[1] && DEBUGGING){
					printf("Merge Error \n");
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
			if(DEBUGGING &&( (isnan(Reps.weights[keptInd[0]])) || (isnan(Reps.weights[keptInd[1]])))){
				printf("WARNING: Moving NAN Weight \n");
			}
			Reps.weights[keptInd[0]] = Reps.weights[keptInd[0]] + Reps.weights[keptInd[1]];

			/*Replace the old simulation with the final non-NAN simulation*/

			if(DEBUGGING && isnan(Reps.weights[Reps.iSimMax])){
				printf("WARNING: Moving NAN Weight \n");
			}

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
			smolFreeSim(Reps.sims[keptInd[1]]);
			Reps.sims[keptInd[1]] = Reps.sims[Reps.iSimMax];
			Reps.weights[keptInd[1]] = Reps.weights[Reps.iSimMax];
			Reps.binLocs[keptInd[1]] = Reps.binLocs[Reps.iSimMax];

			//Setting things to NAN is technically unnecessary. Freeing sim should save memory.
			Reps.sims[Reps.iSimMax] = NULL;
			Reps.weights[Reps.iSimMax] = NAN;
			Reps.binLocs[Reps.iSimMax] = NAN;
			Reps.iSimMax--;

			/*Reorganize the binContents matrix*/
			Reps.binContents[rowCol[2]][mergeBin] = Reps.binContents[Reps.binContentsMax[mergeBin] -1][mergeBin];
			Reps.binContents[-1+Reps.binContentsMax[mergeBin]][mergeBin] = NAN;
			Reps.binContentsMax[mergeBin]--;

			for(entryCheck1 = 0; entryCheck1 < Reps.binContentsMax[mergeBin]-1;entryCheck1++){
				for(entryCheck2 = entryCheck1 + 1; entryCheck2< Reps.binContentsMax[mergeBin]; entryCheck2++){
					if(DEBUGGING && Reps.binContents[entryCheck1][mergeBin]==Reps.binContents[entryCheck2][mergeBin]){
						printf("ERROR: Duplicate entries in BC \n");
					}
				} // finished entryCheck2 loop through this bin's contents
			} // finished entryCheck1 loop through this bin's contents
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
			copySim1(splitInd,Reps.iSimMax+1);
			Reps.weights[splitInd] = Reps.weights[splitInd] / 2;
			Reps.weights[Reps.iSimMax+1] = Reps.weights[splitInd];
			Reps.binLocs[Reps.iSimMax+1] = Reps.binLocs[splitInd];
			Reps.iSimMax++;
			if(DEBUGGING && Reps.iSimMax > ISIMMAXMAX){
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

	if(DEBUGGING && fabs(weightSum-1)>1e-6){
		printf("Weight not conserved \n");
	}

	// loop through replicas in flux bin and delete them
	if(nFlux>0){
		for(iReps = nFlux -1; iReps >=0; iReps--){
				fluxOut += Reps.weights[Reps.binContents[iReps][paramsWe.fluxBin]];
				// Create new replica to replace the one recently lost
				Reps.sims[Reps.binContents[iReps][paramsWe.fluxBin]] = Reps.sims[Reps.iSimMax];
				Reps.weights[Reps.binContents[iReps][paramsWe.fluxBin]] = Reps.weights[Reps.iSimMax];
				Reps.binLocs[Reps.binContents[iReps][paramsWe.fluxBin]] = Reps.binLocs[Reps.iSimMax];

				for(iSimMaxReplace = 0; iSimMaxReplace < Reps.binContentsMax[Reps.binLocs[Reps.iSimMax]];iSimMaxReplace++){
					if(Reps.binContents[iSimMaxReplace][Reps.binLocs[Reps.iSimMax]] == Reps.iSimMax){
						Reps.binContents[iSimMaxReplace][Reps.binLocs[Reps.iSimMax]] = Reps.binContents[iReps][paramsWe.fluxBin];
					}
				}

				Reps.weights[Reps.iSimMax] = NAN;
				Reps.binLocs[Reps.iSimMax] = NAN;
				Reps.iSimMax--;
				Reps.binContentsMax[paramsWe.fluxBin]--;
		}
	}

	if(DEBUGGING && Reps.binContentsMax[paramsWe.fluxBin] != 0){
			printf("Non Zero weight in flux bin \n");
	};

	// Re-weight the remaining replicas
	if(fluxOut != 0){
			for(iSim = 0; iSim < Reps.iSimMax; iSim++){
				Reps.weights[iSim] = Reps.weights[iSim]/(weightSum-fluxOut);
			}
	}

	return fluxOut;
}

void getParams(FILE *DEFile, FILE *WEFile){
	char tmpStr[32];
	int jBin;
	
	fscanf(WEFile,"%s %i", tmpStr, &paramsWe.tau);
	fscanf(WEFile,"%s %i", tmpStr, &paramsWe.repsPerBin);
	fscanf(WEFile,"%s %i", tmpStr, &paramsWe.nInit);
	fscanf(WEFile,"%s %i", tmpStr, &paramsWe.tauMax);
	fscanf(WEFile,"%s %i", tmpStr, &paramsWe.nBins);
	fscanf(WEFile,"%s %i", tmpStr, &paramsWe.fluxBin);
	fscanf(DEFile, "%s %lf", tmpStr, &paramsDe.dt);
	fscanf(DEFile, "%s %lf", tmpStr, &paramsDe.worldLength);
	fscanf(DEFile, "%s %lf", tmpStr, &paramsDe.roiR);
	fscanf(DEFile, "%s %lf", tmpStr, &paramsDe.difC);
	fscanf(DEFile, "%s %i", tmpStr, &paramsDe.nPart);
	paramsWe.fluxBin = 0;
	Reps.iSimMax = paramsWe.nInit -1;
	fclose(DEFile);
	fclose(WEFile);
	Reps.nBins = paramsWe.nBins;
	
	lowBounds[0] = -paramsDe.worldLength/2;
	lowBounds[1] = -paramsDe.worldLength/2;
	highBounds[0] = paramsDe.worldLength/2;
	highBounds[1] = paramsDe.worldLength/2;
	botLeftCornerRect[0] = -paramsDe.worldLength/2;
	botLeftCornerRect[1] = -paramsDe.worldLength/2;
	botLeftCornerRect[2] = paramsDe.worldLength;
	topRightCornerRect[0] = paramsDe.worldLength/2;
	topRightCornerRect[1] = paramsDe.worldLength/2;
	topRightCornerRect[2] = -paramsDe.worldLength;
	roiParams[0] = 0.0;
	roiParams[1] = 0.0;
	roiParams[2] = paramsDe.roiR;
	roiParams[3] = 30;
	insideRoi[0] = 0.0;
	insideRoi[1] = 0.0;
	
	fluxCDF.Tinit = paramsWe.tauMax / 100;
	fluxCDF.nT = paramsWe.tauMax / 100;
	fluxCDF.fluxMax = 0;
	for(jBin = 0; jBin < NFLUXBINS; jBin++){
		fluxCDF.binCounts[jBin] = 0;
	}
	
	printf("Parameters loaded\n");
}

void fluxPDF(double fluxIn){
	int iBinpdf;
	for(iBinpdf = 0; iBinpdf < NFLUXBINS; iBinpdf++){
		if(fluxIn > fluxCDF.binDefs[iBinpdf] && fluxIn < fluxCDF.binDefs[iBinpdf+1]){
			fluxCDF.binCounts[iBinpdf]++;
		}
	}
}

void fluxBin(){
	int iBin;
	double binStep;
	fluxCDF.binDefs[0] = 0;
	binStep = fluxCDF.fluxMax / NFLUXBINS;
	for(iBin = 1; iBin < NFLUXBINS + 1; iBin++){
		fluxCDF.binCounts[iBin] = (double)iBin*binStep;
	}
}

void KSTest(){
	double cdf1, cdf2, ksStat;
	int iBin;
	ksStat = 0;
	for(iBin = 0; iBin < NFLUXBINS; iBin++){
		cdf1 += (double)fluxCDF.oldCounts[iBin]/();
		cdf2 += (double)fluxCDF.binCounts[iBin]/();
		if (fabs(cdf1-cdf2)>ksStat){
			ksStat = fabs(cdf1-cdf2);
			fluxCDF.KSstat = ksStat;
		}
	}
}

int main(int argc, char *argv[]){
	//argv 1: ending simfile, argv2: flux file, argv3: seed / error file, argv4: save / replace rng bit argv5: Execution time file
	//dynamics params: dt, L, R, D, N
	//WE Params: tau, mTarg, tauMax, nBins, ((flux bin))
	int tauQuarter, tauMax, rngBit, iBin, nWE, iSim;
	double fluxAtStep;
	clock_t start[4], stop[4]; //initialDistTime, splitMergeTime, dynamicsTime, totalTime, this also corresponds to the order written in the output file

	//Load simulation / WE parameters from outside files
	FILE *DEFile, *WEFile, *FLFile, *SIMFile, *errFile, *clockFile;
	DEFile = fopen("dynamicsParams.txt","r");
	WEFile = fopen("WEParams.txt","r");
	errFile = fopen(argv[3], "w");
	
	getParams(DEFile, WEFile);

	tauQuarter = paramsWe.tauMax / 4;
	tauMax = paramsWe.tauMax;

	printf("Tau loops + flux vector made \n");

	//Set + record new RNG seed, set up initial distribution
	rngBit = atoi(argv[4]);
	fprintf(errFile, "iseed=%lx\n", RanInitReturnIseed(rngBit));
    fclose(errFile);
	start[0] = clock();
	start[3] = clock();
	initialDist(paramsWe.nInit);
	stop[0] = clock();
	//Set other key initial values in replicas struct
	for(iSim = 0; iSim < paramsWe.nInit; iSim++){
		Reps.weights[iSim] = (double)1/paramsWe.repsPerBin;
		Reps.binLocs[iSim] = findBin(Reps.sims[iSim]);
		Reps.binContents[Reps.binContentsMax[Reps.binLocs[iSim]]][Reps.binLocs[iSim]] = iSim;
		Reps.binContentsMax[Reps.binLocs[iSim]]++;
	}

	printf("Initial Distribution Made \n");
	printf("Initial Bin Location = %i \n", Reps.binLocs[0]);

	//First quarter simulation
	for(nWE = 0; nWE < tauMax; nWE++){
		if(DEBUGGING){
		printf("Tau Step: %i \n", nWE);
		}

		FLFile = fopen(argv[2],"a");
		fluxAtStep = fluxes();
		fprintf(FLFile, "%E \n", fluxAtStep);
		fluxPDF(fluxAtStep);
		fclose(FLFile);
		if(fluxAtStep > fluxCDF.fluxMax && nWE < fluxCDF.Tinit){
			fluxCDF.fluxMax = 3*fluxAtStep;
			fluxBin();
		}
		
		if(nWE == 10){
			start[1] = clock();
		}
		
		splitMerge(nWE);
		
		if(nWE ==10){
			stop[1] = clock();
		}
		
		if(DEBUGGING){
			printf("iSimMax: %i \n", Reps.iSimMax);
		}
		
		for(int iBin = 0; iBin < Reps.nBins; iBin++){
			Reps.binContentsMax[iBin] = 0;
		}
		for(iSim = 0; iSim < Reps.iSimMax; iSim++){
			if(DEBUGGING){
				printf("\n \n \n sim = %i \n \n \n",iSim);
						 }
			if(nWE == 10){
				start[2] = clock();
			}
			dynamicsEngine(Reps.sims[iSim]);
			if(nWE == 10){
				stop[2] = clock();
			}
			Reps.binLocs[iSim] = findBin(Reps.sims[iSim]);
			Reps.binContents[Reps.binContentsMax[Reps.binLocs[iSim]]][Reps.binLocs[iSim]] = iSim;
			Reps.binContentsMax[Reps.binLocs[iSim]]++;
		}
	}
	stop[3] = clock();

	clockFile = fopen(argv[5],"w");
	fprintf(clockFile,"%E \n %E \n %E \n %E \n",(double)(stop[0]-start[0])/CLOCKS_PER_SEC, (double)(stop[1]-start[1])/CLOCKS_PER_SEC, (double)(stop[2]-start[2])/CLOCKS_PER_SEC, (double)(stop[3]-start[3])/CLOCKS_PER_SEC);

	SIMFile = fopen(argv[1], "w");
	for(iBin = 0;iBin <= Reps.iSimMax; iBin++){
		fprintf(SIMFile, "%i, %E \n", findBin(Reps.sims[iBin]), Reps.weights[iBin]);
	}
	fclose(SIMFile);


	return 0;
}
