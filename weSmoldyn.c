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
	FILE *debugFile;

	double p0; /*Doubles giving the merge probabilities*/
	//unsigned int binConPrev[BINCONTENTSMAXMAX][NBINSMAX];
	/*Merging loop*/
	for(mergeBin = binMin; mergeBin <binMax; mergeBin++){
		/*Initialize the merge indices with the first 2 indices stored in appropriate BinContents row*/

		while((Reps.binContentsMax[mergeBin]>paramsWe.repsPerBin)){ // JUN COMMENT: What is the second condition for?
			mergeInd[0] = Reps.binContents[0][mergeBin];
			mergeInd[1] = Reps.binContents[1][mergeBin];
			rowCol[0] = 0;
			rowCol[1] = 1;
			/*Find the locations of the two smallest weights and combine them together*/
			for(repInBin = 0; repInBin < Reps.binContentsMax[mergeBin];repInBin++){ //NUM_Rows needs to change to an appropriate BCM statement

				dummyInd = Reps.binContents[repInBin][mergeBin];
				/*Check the weight associated with dummy Ind with the weight associated with the two merge indices and replace the smaller of the two with this index (assuming this index != mergeInd[0] or mergeInd[1]
				*/
				if(Reps.weights[mergeInd[0]] <= Reps.weights[mergeInd[1]]){
					if(Reps.weights[dummyInd]<Reps.weights[mergeInd[1]] && dummyInd != mergeInd[0]){
						mergeInd[1] = dummyInd;
						rowCol[1] = repInBin;
					}
				}
				else if(Reps.weights[mergeInd[0]]> Reps.weights[mergeInd[1]]){
					if(Reps.weights[dummyInd]<Reps.weights[mergeInd[0]] && dummyInd != mergeInd[1]){
						mergeInd[0] = dummyInd;
						rowCol[0] = repInBin;
					}
				}
			}
			
			/*
			if(mergeInd[0] == mergeInd[1] && DEBUGGING){
					printf("Merge Error \n");
			}*/

			/*Decide which index to keep.*/
			if(Reps.weights[mergeInd[0]]==0 && Reps.weights[mergeInd[1]]==0){
				p0 = 0.5;
			}
			else{
				p0 = Reps.weights[mergeInd[0]] / (Reps.weights[mergeInd[0]]+Reps.weights[mergeInd[1]]);
			}
		
			if(RAND<p0){
				keptInd[0] = mergeInd[0];
				keptInd[1] = mergeInd[1];
				rowCol[2] = rowCol[1];
			}
			else{
				keptInd[0] = mergeInd[1];
				keptInd[1] = mergeInd[0];
				rowCol[2] = rowCol[0];
			}

			/*Update weight of the kept index*/
			if( (isnan(Reps.weights[keptInd[0]])) || (isnan(Reps.weights[keptInd[1]]))){
				debugFile=fopen("Debug.txt","a");
				fprintf(debugFile,"WARNING: Moving NAN Weight. \n keptInd[0] = %i keptWeight[0] = %E \n keptInd[1] = %i keptWeight[1] = %E \n", keptInd[0], Reps.weights[keptInd[0]],keptInd[1], Reps.weights[keptInd[1]]);
				fclose(debugFile);
			}
			Reps.weights[keptInd[0]] += Reps.weights[keptInd[1]];

			/*Replace the old simulation with the final non-NAN simulation*/

			if(DEBUGGING && isnan(Reps.weights[Reps.iSimMax-1])){
				debugFile=fopen("Debug.txt","a");
				fprintf(debugFile,"WARNING: Moving NAN Weight. \n iSimMax = %i weight %E \n", Reps.iSimMax, Reps.weights[Reps.iSimMax-1]);
				fclose(debugFile);
			}

			/*Find iSimMax in binContents and replace it with keptInd[1]*/
			/*For loop: iterates through the row of the bin ISM is in.
			If: Checks to see if the index of the row we're iterating through is ISM. If it is, it replaces it with KI[1], the deleted index.

			The reason we replace ISM with the deleted index is for the incrementing nature of keeping track of the sims.
			We aren't deleting ISM, we're moving it to the spot where the deleted sim was.

			This has problems when the deleted sim is ISM. When that happens, I will take the last element of the table of contents and
			move it to where ISM is in the table.*/

			for(iSimMaxReplace = 0; iSimMaxReplace < Reps.binContentsMax[Reps.binLocs[Reps.iSimMax-1]];iSimMaxReplace++){
				if(Reps.binContents[iSimMaxReplace][Reps.binLocs[Reps.iSimMax-1]] == Reps.iSimMax-1){
					Reps.binContents[iSimMaxReplace][Reps.binLocs[Reps.iSimMax-1]] = keptInd[1];
					simMaxHolder = iSimMaxReplace;
					iSimMaxReplace = Reps.binContentsMax[Reps.binLocs[Reps.iSimMax-1]];
				}
			}
			/*Pretty sure this does the same thing as line 135
			if (Reps.iSimMax-1 ==keptInd[1]){
				Reps.binContents[simMaxHolder][Reps.binLocs[Reps.iSimMax-1]] = Reps.binContents[Reps.binContentsMax[Reps.binLocs[Reps.iSimMax-1]]-1][Reps.binLocs[Reps.iSimMax-1]] //iSimMax gives max numerical, index is 1 less, same for BC Max;
			}*/
			smolFreeSim(Reps.sims[keptInd[1]]);
			Reps.sims[keptInd[1]] = Reps.sims[Reps.iSimMax-1];
			Reps.weights[keptInd[1]] = Reps.weights[Reps.iSimMax-1];
			Reps.binLocs[keptInd[1]] = Reps.binLocs[Reps.iSimMax-1];

			//Setting things to NAN is technically unnecessary. Freeing sim should save memory.
			Reps.sims[Reps.iSimMax-1] = NULL;
			Reps.weights[Reps.iSimMax-1] = NAN;
			Reps.binLocs[Reps.iSimMax-1] = NAN;
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
			copySim1(splitInd,Reps.iSimMax);
			Reps.weights[splitInd] = Reps.weights[splitInd] / 2;
			Reps.weights[Reps.iSimMax] = Reps.weights[splitInd];
			Reps.binLocs[Reps.iSimMax] = Reps.binLocs[splitInd];
			Reps.binContents[-1+Reps.binContentsMax[splitBin]][splitBin] = Reps.iSimMax;
			Reps.binContentsMax[splitBin]++;
			Reps.iSimMax++;
			if(DEBUGGING && Reps.iSimMax > ISIMMAXMAX){
				printf("ERROR: iSimMax out of bounds");
			}
		}
	} // finished splitting loop
	return;
}

double fluxes(){
	double fluxOut = 0;
	double weightSum;
	int jWeight, iReps, simAReplace, iSim, simA;
	FILE *debugFile;
	int nFlux = Reps.binContentsMax[paramsWe.fluxBin]; // Number of replicas in the fluxBin
	
				if(DEBUGGING){
				debugFile = fopen("Debug.txt","a");
				fprintf(debugFile,"%i sims in flux bin \n", nFlux);
				fclose(debugFile);
			}
	
	weightSum = 0;
	for(jWeight = 0; jWeight < Reps.iSimMax; jWeight++){
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
				for(iSim = Reps.iSimMax -1; iSim>=0 ; iSim--){
					if(Reps.binLocs[iSim] != paramsWe.fluxBin){
						simA = iSim;
						iSim = -1;
					}
				}
				//Pretty sure these next few lines should be uncommented out to deal with edge case where iSimMax is in the flux bin
				
				if(Reps.binContents[iReps][paramsWe.fluxBin]== Reps.iSimMax-1){
					simA=Reps.iSimMax-1;
					Reps.sims[simA] = NULL;
					Reps.weights[simA] = NAN;
					Reps.binLocs[simA] = paramsWe.fluxBin;
				}
				else{
				Reps.sims[Reps.binContents[iReps][paramsWe.fluxBin]] = Reps.sims[simA];
				Reps.weights[Reps.binContents[iReps][paramsWe.fluxBin]] = Reps.weights[simA];
				Reps.binLocs[Reps.binContents[iReps][paramsWe.fluxBin]] = findBin(Reps.sims[Reps.binContents[iReps][paramsWe.fluxBin]]);

				for(simAReplace = 0; simAReplace < Reps.binContentsMax[Reps.binLocs[simA]];simAReplace++){
					if(Reps.binContents[simAReplace][Reps.binLocs[simA]] == simA){
						Reps.binContents[simAReplace][Reps.binLocs[simA]] = Reps.binContents[iReps][paramsWe.fluxBin];
					}
				}
				Reps.weights[simA] = NAN;
				Reps.binLocs[simA] = paramsWe.fluxBin;
				Reps.sims[simA] = NULL;
				}

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
	Reps.iSimMax = paramsWe.nInit;
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
	for(jBin = 0; jBin < NFLUXBINS; jBin++){
		fluxCDF.binCounts[jBin] = 0;
	}
	
	printf("Parameters loaded\n");
}

void loadBinDefs(FILE *binParams, FILE *binDefinitions){
	
	char tmpChar[32];
	int iBin, iNbin;
	fscanf(binParams,"%s %i", tmpChar, &binDefs.customBins);
	fscanf(binParams,"%s %i", tmpChar, &binDefs.currentDims);
	fscanf(binParams,"%s %i",tmpChar, &binDefs.nBins);
	if(binDefs.customBins == 1){
	for(int iBinDef = 0; iBinDef < (2*binDefs.currentDims*binDefs.nBins); iBinDef++){
		fscanf(binDefinitions,"%lf", &binDefs.binDefArray[iBinDef]);
	}
	}
}

void KSTest(FILE *FluxFile, int nWE){
	FILE *KSFile; //, *debugKS;
	int iLine, iBin;
	double binStep, cdf1, cdf2, ksStat;
	int nLines = nWE;
	double fluxVect[nWE];
	
	//Clear previous binCounts
	for(iBin = 0; iBin < NFLUXBINS; iBin++){
		fluxCDF.binCounts[iBin] = 0;
		fluxCDF.oldCounts[iBin] = 0;
	}

	//load flux
	fluxCDF.fluxMax = 0;
	for(iLine = 0; iLine < nLines; iLine++){
		fscanf(FluxFile, "%lE\n",&fluxVect[iLine]);
		if(iLine> nLines/3 && fluxVect[iLine]>fluxCDF.fluxMax){
			fluxCDF.fluxMax = fluxVect[iLine];
		}
	}

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
	fclose(KSFile);
	fluxCDF.nT *= 2;
	return;
}


int main(int argc, char *argv[]){
	//argv 1: ending simfile, argv2: flux file, argv3: seed / error file, argv4: save / replace rng bit argv5: Execution time file
	//dynamics params: dt, L, R, D, N
	//WE Params: tau, mTarg, tauMax, nBins, ((flux bin))
	int tauMax, rngBit, iBin, nWE, iSim, iBCM, nanCheck, firstNAN, iDimer,iBinContents; //tauQuarter omitted
	double fluxAtStep, binWeight;
	clock_t start[4], stop[4]; //initialDistTime, splitMergeTime, dynamicsTime, totalTime, this also corresponds to the order written in the output file

	//Load simulation / WE parameters from outside files
	FILE *DEFile, *WEFile, *FLFile, *SIMFile, *errFile, *clockFile, *mCountsFile, *structStoreFile, *structNANFile, *debugFile, *mCountsWeightedFile;
	char *fluxFileStr;
	FILE *binDefinitions, *binParams;
	fluxFileStr = argv[2];
	DEFile = fopen("dynamicsParams.txt","r");
	WEFile = fopen("WEParams.txt","r");
	errFile = fopen(argv[3], "w");
	
	getParams(DEFile, WEFile);
	fclose(DEFile);
	fclose(WEFile);
	binParams = fopen("binParams.txt","r");
	binDefinitions = fopen("binDefinitions.txt","r");
	loadBinDefs(binParams, binDefinitions);
	fclose(binParams);
	fclose(binDefinitions);
	
	//Defining variables for continuous dimer counting, gives more data points for comparing monomerization fraction
	int mCounts[3];
	double  mCountsWeighted[3]; 
	//These files will each give the total monomers, dimers, and either all "time" elapsed for all sims or all weight measured across
	for(iDimer = 0; iDimer < 3; iDimer++){
		mCounts[iDimer] = 0;
		mCountsWeighted[iDimer] = 0;
	}
	
	tauMax = paramsWe.tauMax;
	if(DEBUGGING){
		debugFile = fopen("Debug.txt","a");
		fprintf(debugFile,"Tau loops + flux vector made \n");
		fclose(debugFile);
	}
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
	if(DEBUGGING){
	debugFile = fopen("Debug.txt","a");
	fprintf(debugFile, "Initial Distribution Made \n");
	fclose(debugFile);
	}
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

	start[3] = clock();
	
	//Simulation Loop
	for(nWE = 0; nWE < tauMax; nWE++){
		
		//Print current tau step to stdout for interactive debugging
		if(DEBUGGING){
		debugFile = fopen("Debug.txt","a");
		fprintf(debugFile,"Tau Step: %i \n", nWE);
		fclose(debugFile);
		}
		//Flux Recording
		if(paramsWe.fluxBin >= 0){
		FLFile = fopen(fluxFileStr,"a");
		fluxAtStep = fluxes();
		fprintf(FLFile, "%E\n", fluxAtStep);
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
		if(nWE == fluxCDF.nT&&paramsWe.fluxBin >= 0){
			if(DEBUGGING){
				debugFile = fopen("Debug.txt","a");
				fprintf(debugFile,"KS Recording \n");
				fclose(debugFile);
			}
			FLFile = fopen(fluxFileStr,"r");
			KSTest(FLFile, nWE);
			fclose(FLFile);
			if(DEBUGGING){
				debugFile = fopen("Debug.txt","a");
				fprintf(debugFile,"KS Recording \n");
				fclose(debugFile);
			}
			}
		
		//Storing replicas struct without any stray NANs
		if(DEBUGGING){
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
		if(DEBUGGING){
			if(nanCheck == 0){
			for(iSim = 0; iSim < Reps.iSimMax; iSim++){
				if(Reps.weights[iSim] != Reps.weights[iSim]){
					nanCheck = 1;
				}
			}
			}
		
		
				debugFile = fopen("Debug.txt","a");
				fprintf(debugFile,"Nans checked \n");
				fclose(debugFile);
		
		
		//Recording stray NANs (only the first time they appear)
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
		}
		
		if(nWE ==10){
			stop[1] = clock();
		}
		
		if(DEBUGGING){
			debugFile = fopen("Debug.txt","a");
			fprintf(debugFile,"NANs checked and recorded. iSimMax = %i \n Dynamics Starting \n", Reps.iSimMax);
			fclose(debugFile);
		}
		
		//Dynamics and rebuilding BCM
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
			dynamicsEngine(Reps.sims[iSim], mCounts, mCountsWeighted, nWE, iSim);
			
			if(nWE == 10){
				stop[2] = clock();
			}
			Reps.binLocs[iSim] = findBin(Reps.sims[iSim]);
			Reps.binContents[Reps.binContentsMax[Reps.binLocs[iSim]]][Reps.binLocs[iSim]] = iSim;
			Reps.binContentsMax[Reps.binLocs[iSim]]++;
			//Dimerization Fraction Recording
			if(!MONOFRACEACHDT){
 			if(nWE > tauMax/2){
				mCounts[0] += Reps.sims[iSim]->mols->nl[0];
				mCounts[1] += Reps.sims[iSim]->mols->nl[1];
				mCounts[2]++;
				mCountsWeighted[0] += (double) Reps.sims[iSim]->mols->nl[0] * Reps.weights[iSim];
				mCountsWeighted[1] += (double) Reps.sims[iSim]->mols->nl[1] * Reps.weights[iSim];
				mCountsWeighted[2] += Reps.weights[iSim];
			} 
			}
	}
		
		//Record weight distribution among bins
		/*
		if(tauMax >= 1000){
		if( (nWE+1) % (tauMax/1000) == 0 || nWE < 1000){
		SIMFile = fopen(argv[1],"a");
		for(iBin = 0; iBin < Reps.nBins; iBin++){
			binWeight = 0;
			for(iBinContents=0;iBinContents < Reps.binContentsMax[iBin];iBinContents++){
				binWeight += Reps.weights[Reps.binContents[iBinContents][iBin]];
			}
			fprintf(SIMFile,"%i, %E \n",iBin, binWeight);
		}
		fprintf(SIMFile, "\n"); //Consider changing this. Can also consider making new file for each tau step
		//overall sizes should stay reasonable either way
		fclose(SIMFile);
		
		}}
		else{
		SIMFile = fopen(argv[1],"a");
		for(iBin = 0; iBin < Reps.nBins; iBin++){
			binWeight = 0;
			for(iBinContents=0;iBinContents < Reps.binContentsMax[iBin];iBinContents++){
				binWeight += Reps.weights[Reps.binContents[iBinContents][iBin]];
			}
			fprintf(SIMFile,"%i, %E \n",iBin, binWeight);
		}
		fprintf(SIMFile, "\n"); //Consider changing this. Can also consider making new file for each tau step
		//overall sizes should stay reasonable either way
		fclose(SIMFile);			
		}*/
		
		if(DEBUGGING){
			debugFile = fopen("Debug.txt","a");
			fprintf(debugFile,"Dynamics Finished  \n Recording mCounts \n");
			fclose(debugFile);
		}
		
		if(nWE > tauMax/2){
			//mCountsFile = fopen("mCounts.txt","a");
			//fprintf(mCountsFile, "%i, %i, %i \n", mCounts[0], mCounts[1], mCounts[2]);
			//fclose(mCountsFile);
			mCountsWeightedFile = fopen("mCountsWeighted.txt","a");
			fprintf(mCountsWeightedFile, "%f, %f, %f \n", mCountsWeighted[0], mCountsWeighted[1],mCountsWeighted[2]);
			fclose(mCountsWeightedFile);
		}
	}
	
	stop[3] = clock();

	
	//Time Recording
	clockFile = fopen(argv[5],"w");
	fprintf(clockFile,"%E \n %E \n %E \n %E \n",(double)(stop[0]-start[0])/CLOCKS_PER_SEC, (double)(stop[1]-start[1])/CLOCKS_PER_SEC, (double)(stop[2]-start[2])/CLOCKS_PER_SEC, (double)(stop[3]-start[3])/CLOCKS_PER_SEC);
	fclose(clockFile);
	
	//Free Memory and finish
	freeAllSims();
	return 0;
}