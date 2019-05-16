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
	int splitBin, mergeBin, repInBin, iSimMaxReplace; //Looping variables. entryCheck1 and 2 should be declared here if checking

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

			/*for(entryCheck1 = 0; entryCheck1 < Reps.binContentsMax[mergeBin]-1;entryCheck1++){
				for(entryCheck2 = entryCheck1 + 1; entryCheck2< Reps.binContentsMax[mergeBin]; entryCheck2++){
					if(DEBUGGING && Reps.binContents[entryCheck1][mergeBin]==Reps.binContents[entryCheck2][mergeBin]){
						printf("ERROR: Duplicate entries in BC \n");
					}
				} // finished entryCheck2 loop through this bin's contents
			} // finished entryCheck1 loop through this bin's contents */
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
	int jWeight, iReps, simAReplace, iSim, simA;
	int nFlux = Reps.binContentsMax[paramsWe.fluxBin]; // Number of replicas in the fluxBin

	weightSum = 0;
	for(jWeight = 0; jWeight <= Reps.iSimMax; jWeight++){
		weightSum = weightSum + Reps.weights[jWeight];
	}

	if(DEBUGGING && fabs(weightSum-1)>1e-6){
		printf("Weight not conserved \n");
	}

	// loop through replicas in flux bin and delete them. Replace indices with empty non-flux indices
	if(nFlux>0){
		for(iReps = nFlux -1; iReps >=0; iReps--){
				fluxOut += Reps.weights[Reps.binContents[iReps][paramsWe.fluxBin]];
			
				//Delete simulation. scan through all sims, find the last index that isn't in flux bin, then call that simA
				//Move simA pointer to weights to iReps location 
				//free the sim of deleted sim too
				smolFreeSim(Reps.sims[Reps.binContents[iReps][paramsWe.fluxBin]]);
				for(iSim = Reps.iSimMax; iSim >= 0; iSim--){
					if(Reps.binLocs[iSim] != paramsWe.fluxBin){
						simA = iSim;
					}
				}
				Reps.sims[Reps.binContents[iReps][paramsWe.fluxBin]] = Reps.sims[simA];
				Reps.weights[Reps.binContents[iReps][paramsWe.fluxBin]] = Reps.weights[simA];
				Reps.binLocs[Reps.binContents[iReps][paramsWe.fluxBin]] = Reps.binLocs[simA];

				for(simAReplace = 0; simAReplace < Reps.binContentsMax[Reps.binLocs[simA]];simAReplace++){
					if(Reps.binContents[simAReplace][Reps.binLocs[simA]] == simA){
						Reps.binContents[simAReplace][Reps.binLocs[simA]] = Reps.binContents[iReps][paramsWe.fluxBin];
					}
				}
				Reps.weights[simA] = NAN;
				Reps.binLocs[simA] = paramsWe.fluxBin;
		}
	}
	for(iReps = 0; iReps < nFlux; iReps++){
				Reps.iSimMax--;
				Reps.binContentsMax[paramsWe.fluxBin]--;
	}
	
	if( Reps.binContentsMax[paramsWe.fluxBin] != 0){
		if(DEBUGGING){
			printf("BCM != 0 after flux step\n");
		}
		Reps.binContentsMax[paramsWe.fluxBin] = 0;
	}
	
		
	//

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
	fscanf(DEFile, "%s %lf", tmpStr, &paramsDe.difM);
	fscanf(DEFile, "%s %lf", tmpStr, &paramsDe.difD);
	fscanf(DEFile, "%s %lf", tmpStr, &paramsDe.bindR);
	fscanf(DEFile, "%s %lf", tmpStr, &paramsDe.unbindK);
	fscanf(DEFile, "%s %i", tmpStr, &paramsDe.nPart);
	fscanf(DEFile, "%s %i", tmpStr, &paramsDe.reactBit);
	//paramsWe.fluxBin = 0;
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
	
	fluxCDF.nT = 300;
	fluxCDF.fluxMax = 0.001;
	for(jBin = 0; jBin < NFLUXBINS; jBin++){
		fluxCDF.binCounts[jBin] = 0;
	}
	
	printf("Parameters loaded\n");
}

void KSTest(char *flFile){
	FILE *KSFile, *FluxFile;
	char chrScan;
	int iLine, iBin;
	double binStep, cdf1, cdf2, ksStat;
	int nLines = 0;
	
	//Clear previous binCounts
	for(iBin = 0; iBin < NFLUXBINS; iBin++){
		fluxCDF.binCounts[iBin] = 0;
		fluxCDF.oldCounts[iBin] = 0;
	}
	
	//open files for reading
	FluxFile = fopen(flFile,"r");
	
	//count number of flux measurements
	chrScan = getc(FluxFile);
	while(chrScan != EOF){
		if(chrScan == 'n'){
			nLines++;
		}
		chrScan = getc(FluxFile);
	}
	//declare variable and load flux into it
	double fluxVect[nLines];
	fluxCDF.fluxMax = 0;
	for(iLine = (nLines/3)-1; iLine < nLines; iLine++){
		fscanf(FluxFile, "%E",fluxVect[iLine]);
		if(fluxVect[iLine]>fluxCDF.fluxMax){
			fluxCDF.fluxMax = fluxVect[iLine];
		}
	}
	fclose(FluxFile);
	
	
	//Create new Bin Defs
	binStep = fluxCDF.fluxMax / NFLUXBINS;
	fluxCDF.binDefs[0] = 0;
	for(iBin = 1; iBin <= NFLUXBINS; iBin++){
		fluxCDF.binDefs[iBin] = (double)iBin*binStep;
	}
	
	//Update Bin Counts
	for(iLine = (nLines/3) - 1; iLine < (2*nLines/3); iLine++){
		for(iBin=0;iBin<NFLUXBINS;iBin++){
			if(fluxVect[iLine]>fluxCDF.binDefs[iBin] && fluxVect[iLine]<fluxCDF.binDefs[iBin+1]){
				fluxCDF.oldCounts[iBin]++;
			}
		}
	}
	
	for(iLine = (2*nLines/3); iLine < nLines; iLine++){
		for(iBin=0;iBin<NFLUXBINS;iBin++){
			if(fluxVect[iLine]>fluxCDF.binDefs[iBin] && fluxVect[iLine]<fluxCDF.binDefs[iBin+1]){
				fluxCDF.binCounts[iBin]++;
			}
		}
	}
	
	//Perform KS Test
	cdf1 = 0;
	cdf2 = 0;
	ksStat = 0;
	for(iBin = 0; iBin < NFLUXBINS; iBin++){
		cdf1 += (double)fluxCDF.oldCounts[iBin]/(fluxCDF.nT/3);
		cdf2 += (double)fluxCDF.binCounts[iBin]/(fluxCDF.nT/3);
		if(fabs(cdf1-cdf2)>ksStat){
			ksStat = fabs(cdf1-cdf2);
		}
	}
	
	KSFile = fopen("ksOut.txt", "a");
	fprintf(KSFile, "maxFlux = %E \n", fluxCDF.fluxMax);
	for(iBin = 0; iBin < NFLUXBINS; iBin++){
		fprintf(KSFile, "%i ", fluxCDF.binCounts[iBin]);
	}
	fprintf(KSFile, "\n");
	for(iBin = 0; iBin < NFLUXBINS; iBin++){
		fprintf(KSFile, "%i ", fluxCDF.oldCounts[iBin]);
	}
	fprintf(KSFile, "\n ksStat: %E \n nT %i \n",ksStat, fluxCDF.nT);
	fluxCDF.nT *= 2;
}


int main(int argc, char *argv[]){
	//argv 1: ending simfile, argv2: flux file, argv3: seed / error file, argv4: save / replace rng bit argv5: Execution time file
	//dynamics params: dt, L, R, D, N
	//WE Params: tau, mTarg, tauMax, nBins, ((flux bin))
	int tauMax, rngBit, iBin, nWE, iSim, iBCM, nanCheck, firstNAN, iDimer; //tauQuarter omitted
	double fluxAtStep;
	clock_t start[4], stop[4]; //initialDistTime, splitMergeTime, dynamicsTime, totalTime, this also corresponds to the order written in the output file

	//Load simulation / WE parameters from outside files
	FILE *DEFile, *WEFile, *FLFile, *SIMFile, *errFile, *clockFile, *ksFile, *mCountsFile, *structStoreFile, *structNANFile, *debugFile;
	char *fluxFileStr;
	fluxFileStr = argv[2];
	DEFile = fopen("dynamicsParams.txt","r");
	WEFile = fopen("WEParams.txt","r");
	errFile = fopen(argv[3], "w");
	
	getParams(DEFile, WEFile);
	
	//Defining double for continuous dimer counting, gives more data points for comparing monomerization fraction
	int dMax = paramsDe.nPart/2;
	double dCounts[dMax + 2]; //Max dimers = nPart/2, min dimers = 0; so nPart/2 +1 possibilities for dimer frac, put in 1 more for storing number of recordings
	for(iDimer = 0; iDimer < dMax +2; iDimer++){
		dCounts[iDimer] = 0;
	}
	
	tauMax = paramsWe.tauMax;
	debugFile = fopen("Debug.txt","a");
	fprintf(debugFile,"Tau loops + flux vector made \n");
	fclose(debugFile);
	//Set + record new RNG seed, set up initial distribution
	rngBit = atoi(argv[4]);
	fprintf(errFile, "iseed=%lx\n", RanInitReturnIseed(rngBit));
	fprintf(errFile, "fluxBin = %i \n", paramsWe.fluxBin);
	if(!STOPCOMMAND){
		fprintf(errFile, "Simulations do NOT stop upon evacuation \n");
	}
	if(!WEENABLE){
		fprintf(errFile, "Weighted Ensemble process is deactivated \n");
	}
	if(paramsWe.fluxBin < 0){
		fprintf(errFile, "Flux measurement disabled \n");
	}
    fclose(errFile);
	start[0] = clock();
	initialDist(paramsWe.nInit);
	stop[0] = clock();
	debugFile = fopen("Debug.txt","a");
	fprintf(debugFile, "Initial Distribution Made \n");
	fclose(debugFile);
	//Set other key initial values in replicas struct
	for(int iBin = 0; iBin < Reps.nBins; iBin++){
			Reps.binContentsMax[iBin] = 0;
	}
	for(iSim = 0; iSim < paramsWe.nInit; iSim++){
		Reps.weights[iSim] = (double)1/paramsWe.nInit;
		Reps.binLocs[iSim] = findBin(Reps.sims[iSim]);
		Reps.binContents[Reps.binContentsMax[Reps.binLocs[iSim]]][Reps.binLocs[iSim]] = iSim;
		Reps.binContentsMax[Reps.binLocs[iSim]]++;
	}
	nanCheck = 0;
	firstNAN = 0;

	printf("Initial Distribution Made \n");
	printf("Initial Bin Location = %i \n", Reps.binLocs[0]);

	start[3] = clock();
	
	//Simulation Loop
	for(nWE = 0; nWE < tauMax; nWE++){
		
		//Print current tau step to stdout for interactive debugging
		debugFile = fopen("Debug.txt","a");
		fprintf(debugFile,"Tau Step: %i \n", nWE);
		fclose(debugFile);
		
		//Flux Recording
		if(paramsWe.fluxBin >= 0){
		FLFile = fopen(fluxFileStr,"a");
		fluxAtStep = fluxes();
		fprintf(FLFile, "%E \n", fluxAtStep);
		fclose(FLFile);
			if(DEBUGGING){
				debugFile = fopen("Debug.txt","a");
				fprintf(debugFile,"Flux Recorded \n");
				fclose(debugFile);
			}
		}
		
		//Single Step Time Measurement
		if(nWE == 10){
			start[1] = clock();
		}
		
		//KS Recording
		if(nWE == fluxCDF.nT){
			KSTest(fluxFileStr);
			}
		
		//Storing replicas struct without any stray NANs
		if(nanCheck == 0){
			structStoreFile = fopen("SimStructs.txt", "w");
			fprintf(structStoreFile, "Tau = %i \n", nWE);
			for(iSim = 0; iSim < Reps.iSimMax; iSim++){
				fprintf(structStoreFile, "%i %E \n", Reps.binLocs[iSim], Reps.weights[iSim]);
			}
			for(iBin = 0; iBin < NBINSMAX; iBin++){
				for(iBCM = 0; iBCM < Reps.binContentsMax[iBin]; iBCM++){
					fprintf(structStoreFile, "%i ", Reps.binContents[iBCM][iBin]);
				}
				fprintf(structStoreFile, "\n");
			}
			fclose(structStoreFile);
		}
		
		//Splitting and Merging
		if(WEENABLE){
			if(DEBUGGING){
				debugFile = fopen("Debug.txt","a");
				fprintf(debugFile,"Split Merge Starting... \n");
				fclose(debugFile);
			}
		splitMerge(nWE);
			if(DEBUGGING){
				debugFile = fopen("Debug.txt","a");
				fprintf(debugFile,"Split Merge completed \n Checking for NANs \n");
				fclose(debugFile);
			}
		}
		
		
		//Checking for first appearance of stray NANs
		/*for(iSim = 0; (iSim < Reps.iSimMax || nanCheck == 0); iSim++){
			if(Reps.weights[iSim] != Reps.weights[iSim]){
				nanCheck = 1;
			}
		}*/
		
		if(DEBUGGING){
				debugFile = fopen("Debug.txt","a");
				fprintf(debugFile,"Nans checked \n");
				fclose(debugFile);
			}
		
		//Recording stray NANs (only the first time they appear)
		/*if(nWE == 0){
		if((nanCheck == 1 && firstNAN == 0)){
			structNANFile = fopen("nanStructs.txt","w");
			for(iSim = 0; iSim < Reps.iSimMax; iSim++){
				fprintf(structNANFile, "%i %E \n", Reps.binLocs[iSim],Reps.weights[iSim]);
			}
			for(iBin = 0; iBin < NBINSMAX; iBin++){
				for(iBCM = 0; iBCM < Reps.binContentsMax[iBin]; iBCM++){
					fprintf(structNANFile, "%i ", Reps.binContents[iBCM][iBin]);
				}
				fprintf(structNANFile, "\n");
			}
			fclose(structNANFile);
			firstNAN = 1;
		}
		}*/
		
		if(nWE ==10){
			stop[1] = clock();
		}
		
		if(DEBUGGING){
			debugFile = fopen("Debug.txt","a");
			fprintf(debugFile,"NANs checked and recorded. iSimMax = %i \n Dynamics Starting \n", Reps.iSimMax);
			fclose(debugFile);
		}
		
		for(int iBin = 0; iBin < Reps.nBins; iBin++){
			Reps.binContentsMax[iBin] = 0;
		}
		for(iSim = 0; iSim < Reps.iSimMax; iSim++){
			/*if(DEBUGGING){
				debugFile = fopen("Debug.txt","a");
				fprintf(debugFile,"Sim %i dynamics starting  \n", iSim);
				fclose(debugFile);
						 }*/
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
			dCounts[Reps.sims[iBin]->mols->nl[1]] += Reps.weights[iSim];
		}
		dCounts[dMax+1]++;
		debugFile = fopen("Debug.txt","a");
		fprintf(debugFile,"Dynamics Finished  \n");
		fclose(debugFile);
	}
	
	stop[3] = clock();

	
	//Time Recording
	clockFile = fopen(argv[5],"w");
	fprintf(clockFile,"%E \n %E \n %E \n %E \n",(double)(stop[0]-start[0])/CLOCKS_PER_SEC, (double)(stop[1]-start[1])/CLOCKS_PER_SEC, (double)(stop[2]-start[2])/CLOCKS_PER_SEC, (double)(stop[3]-start[3])/CLOCKS_PER_SEC);
	fclose(clockFile);

	//Final Equilibrium Recording
	debugFile = fopen("Debug.txt","a");
	fprintf(debugFile,"Recording final sims  \n");
	fclose(debugFile);
	SIMFile = fopen(argv[1], "w");
	for(iBin = 0;iBin <= Reps.iSimMax; iBin++){
		fprintf(SIMFile, "%i, %E \n", findBin(Reps.sims[iBin]), Reps.weights[iBin]);
	}
	fclose(SIMFile);

	//Molecule counts recording
	debugFile = fopen("Debug.txt","a");
	fprintf(debugFile,"Recording mCounts  \n");
	fclose(debugFile);
	
	
	
	mCountsFile = fopen("mCounts.txt", "a");
	for(iDimer = 0; iDimer <= dMax; iDimer++){
		fprintf(mCountsFile, "%i, %E \n", iDimer, dCounts[iDimer]); //Change later to not have hard code
	}
	fprintf(mCountsFile, "%E, %E", dCounts[dMax+1], (double)paramsDe.dt*paramsWe.tau);
	fclose(mCountsFile);
	
	//Free Memory and finish
	freeAllSims();
	return 0;
}