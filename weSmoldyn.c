#include "weSmoldyn.h"
#include "smolDynamics.c"

void splitMerge(int nWE){
/*
!---------------------------------------------------------------------------------------------------------------------
!		Description:
!			This function performs WE splitting and merging
!		Notes:
!			iSimMax > 0, so when being used as an index we deduct 1 from its value.
!		Method:
!			1. Define variables for tracking indices. dummyInd + rowCol + mergeInd + keptInd + splitInd. Fairly useful
!				for avoiding very extended index terms (e.g. Reps.Sims[Reps.binContents[Reps.iSimMax][mergeBin]]) as
!				well as some looping variables
!			2. Begin a loop through all bins to check for merging. While a bin has more replicas inside than the mTarg
!				do the following:
!				a. Loop through the replicas in the bin to find the two replicas with the smallest WE weights.
!				b. Determine which index should be kept probabilistically, determined by relative weights. Update
!					dummy indices.
!				c. Combine the weights from both indices into the kept index
!				d. Replace deleted index with another index and update the table of contents matrix (binContents)
!					i. Find iSimMax in binContents, replace its index with the index of the deleted replica
!					ii. Free the sim of the deleted replica, replace the pointer+weight+bin with iSimMax's
!					iii. Place NAN / NULL values where iSimMax's old location was, decrease iSimMax by 1
!					iv. In the BC matrix for the merging bin, replace the deleted index with the index at the end
!						take the index at the end and set to NAN, decrease BC max for the bin.
!			3. Begin a loop through all bins to check for splitting. While a bin has less replicas inside than mTarg
!				but also > 0 replicas, do the following:
!				a. Loop through all replicas in the bin to find the replica with the greateset WE weight.
!				b. Execute copySim1, copying the greatest weight to the first empty simptr in Reps.sims
!				c. Halve the weight of the split sim, copy its weight + bin location to appropriate Reps.weights and 
!					Reps.binLocs
!				d. Add the new index to the BC matrix, increment BC max + iSimMax to account for extra replica
!---------------------------------------------------------------------------------------------------------------------
*/
	int binMin = 0;
	int binMax = Reps.nBins;
	unsigned int dummyInd;/*Dummy index that will hold an index location in both split and merge loop*/
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
	for(mergeBin = binMin; mergeBin <=binMax; mergeBin++){
		/*Initialize the merge indices with the first 2 indices stored in appropriate BinContents row*/

		while((Reps.binContentsMax[mergeBin]>paramsWe.repsPerBin)){
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
				//Occaisionally, for bins with very low total weight, multiple replicas will have weight 0.
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
/*
!---------------------------------------------------------------------------------------------------------------------
!		Description:
!			Measures weight inside flux bin, frees sims in flux bin, readjusts weights of all other sims to sum to 1.
!			The method of reintroducing the flux weight is dependent on the parameter paramsWe.replaceFluxSims
!		Method:
!			Counts number of replicas in flux bin. If there is at least one, do the following:
!				1. First measures total weight in the system and if the DEBUGGING flag is true then it will output 
!					when the weight inside the system deviates from 1 by more than 1e-6.
!				2. If there are replicas in the flux bin, begin the weight replacement process.
!					a. If fluxSims are replaced, then each replica that gets replaced has it's flux added to the total
!						for the step. The replica is then freed and replaced with a freshly initialized replica.
!					b. If flux sims are not replaced, the sims are freed and the weight of the remaining replicas is
!						scaled by the flux removed.
!				3. The binContents matrix and other indexing is tracked in the if/else statements inside the
!					if(nFlux>0) loop. As the number of replicas, Reps.iSimMax, changes in this step by a varying
!					amount, reorganizing the binContents matrix and binLocs matrices requires care. 
!						E.g. Much of the code uses iSimMax as an index of a specific replica. If iSimMax = n
!						if iSimMax-10 is in the	flux bin and removed, iSimMax reduces to n-1
!						and now the replica at n needs to be moved to an index j<=n-1
!---------------------------------------------------------------------------------------------------------------------
*/
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
		if(paramsWe.replaceFluxSims==1){
		for(iReps = 0; iReps < Reps.binContentsMax[paramsWe.fluxBin];iReps++){
			fluxOut += Reps.weights[Reps.binContents[iReps][paramsWe.fluxBin]];
			smolFreeSim(Reps.sims[Reps.binContents[iReps][paramsWe.fluxBin]]);
			buildSingleSim(Reps.binContents[iReps][paramsWe.fluxBin]);
		}
		}
		else{
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
				//Deals with edge case where iSimMax is in the flux bin
				
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
	}
	if(paramsWe.replaceFluxSims==0){
	for(iReps = 0; iReps < nFlux; iReps++){
		
				Reps.iSimMax--;
				Reps.binContentsMax[paramsWe.fluxBin]--;
		}
	}
	
	if( Reps.binContentsMax[paramsWe.fluxBin] != 0){
		if(DEBUGGING){
			printf("BCM != 0 after flux step\n");
		}
		Reps.binContentsMax[paramsWe.fluxBin] = 0;
	}

	// Re-weight the remaining replicas
	if(fluxOut != 0 && paramsWe.replaceFluxSims != 1){
			for(iSim = 0; iSim < Reps.iSimMax; iSim++){
				Reps.weights[iSim] = Reps.weights[iSim]/(weightSum-fluxOut);
			}
	}

	return fluxOut;
}
void getParams(FILE *DEFile, FILE *WEFile, FILE *CorralsFile){
	/*
!---------------------------------------------------------------------------------------------------------------------
!		Description:
!			Loads simulation parameters from text files: "dynamicsParams.txt", "WEParams.txt", "corralsParams.txt".
!		Notes:
!			Bin parameters are not handled in these params. Bin parameters are also separate from WEParams.txt. Bin
!			parameters handled in loadBinDefs(), see below. Many roi parameters are hard coded in, such as the
!			the center of the ROI being at the origin and the origin signifying the region inside the compartment.
!			the value hard-coded as 30 is for drawing purposes. There is no drawing for WE, so this is unimportant.
!		Method:
!			1. Scan text files and give values to global variables defined in weSmoldyn.h.
!---------------------------------------------------------------------------------------------------------------------
*/
	char tmpStr[32];
	int jBin;
	
	fscanf(WEFile,"%s %i", tmpStr, &paramsWe.tau);
	fscanf(WEFile,"%s %i", tmpStr, &paramsWe.repsPerBin);
	fscanf(WEFile,"%s %i", tmpStr, &paramsWe.nInit);
	fscanf(WEFile,"%s %i", tmpStr, &paramsWe.tauMax);
	fscanf(WEFile,"%s %i", tmpStr, &paramsWe.nBins);
	fscanf(WEFile,"%s %i", tmpStr, &paramsWe.fluxBin);
	fscanf(WEFile,"%s %i" ,tmpStr, &fluxCDF.nT);
	fscanf(WEFile,"%s %i", tmpStr, &paramsWe.replaceFluxSims);
	fscanf(DEFile, "%s %lf", tmpStr, &paramsDe.dt);
	fscanf(DEFile, "%s %lf", tmpStr, &paramsDe.worldLength);
	fscanf(DEFile, "%s %lf", tmpStr, &paramsDe.roiR);
	fscanf(DEFile, "%s %lf", tmpStr, &paramsDe.difM);
	fscanf(DEFile, "%s %lf", tmpStr, &paramsDe.difD);
	fscanf(DEFile, "%s %lf", tmpStr, &paramsDe.bindR);
	fscanf(DEFile, "%s %lf", tmpStr, &paramsDe.unbindK);
	fscanf(DEFile, "%s %i", tmpStr, &paramsDe.nPart);
	fscanf(DEFile, "%s %i", tmpStr, &paramsDe.reactBit);
	fscanf(DEFile, "%s %i",tmpStr,&paramsDe.reentryRateBit);
	fscanf(DEFile, "%s %lf",tmpStr, &paramsDe.reentryRate);
	fscanf(DEFile, "%s %i",tmpStr, &paramsDe.densityBit);
	fscanf(DEFile, "%s %lf",tmpStr, &paramsDe.density);
	fscanf(DEFile, "%s %i", tmpStr, &paramsDe.monomerStart);
	fscanf(CorralsFile,"%s %i",tmpStr, &paramsDe.corralsBit);
	fscanf(CorralsFile,"%s %lf",tmpStr, &paramsDe.corralWidth);
	fscanf(CorralsFile,"%s %lf",tmpStr, &paramsDe.corralRate);
	
	if(paramsDe.densityBit){
		paramsDe.worldLength = sqrt((double)paramsDe.nPart/paramsDe.density);
	}
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
	//
	roiParams[0] = 0.0;
	roiParams[1] = 0.0;
	roiParams[2] = paramsDe.roiR;
	roiParams[3] = 30;
	insideRoi[0] = 0.0;
	insideRoi[1] = 0.0;
	
	for(jBin = 0; jBin < NFLUXBINS; jBin++){
		fluxCDF.binCounts[jBin] = 0;
	}
	
	printf("Parameters loaded\n");
}

void loadBinDefs(){
	/*
!---------------------------------------------------------------------------------------------------------------------
!		Description:
!			Uses the "binParams.txt" and "binDefinitions.txt" files to load user-defined bin definitions. 
!		Notes:
!			If there are no custom bins, then the bins will either range from 0 to paramsDe.nPart.
!		Method:
!			1. Load the files
!			2. Change the binDefs struct to account for the custom bin definitions. Update paramsWe.nBins / Reps.nBins
!				to account for this.
!---------------------------------------------------------------------------------------------------------------------
*/
	char tmpChar[32];
	int iBinDef;
	FILE *binDefinitions, *binParams;
	binParams = fopen("binParams.txt","r");
	binDefinitions = fopen("binDefinitions.txt","r");
	fscanf(binParams,"%s %i", tmpChar, &binDefs.customBins);
	fscanf(binParams,"%s %i", tmpChar, &binDefs.currentDims);
	fscanf(binParams,"%s %i",tmpChar, &binDefs.nBins);
	if(binDefs.customBins == 1){
	for(iBinDef = 0; iBinDef < (binDefs.nBins+1); iBinDef++){
		fscanf(binDefinitions,"%lf", &binDefs.binDefArray[iBinDef]);
	}
		paramsWe.nBins = binDefs.nBins;
		Reps.nBins = paramsWe.nBins;
	}
	fclose(binParams);
	fclose(binDefinitions);
}



void KSTest(FILE *FluxFile, int nWE){
	/*
!---------------------------------------------------------------------------------------------------------------------
!		Description:
!			Performs a simplified KS test as well as a modified KS test on the middle third and final third of the data 
!			that has been measured so far. If certain accuracy standard have been met then this will stop the WE sim
!		Notes:
!			Modified KS test involves performing a KS test with "duplicated" zeros removed from both the middle third
!			and the final third of data. The requirement for this KS test is much less strict than the simple KS.
!			If there are not enough non-zero measurements then the results of the KS test will be ignored / the 
!			test will fail automatically.
!		Method:
!			1. Clear old data from previous run of the function
!			2. Load flux measurements from the FluxFile.
!			3. Create histogram bins based on the range of flux values loaded
!			4. Sort flux measurements into histogram bins
!			5. Create modified histogram with excess zeroes removed.
!			6. Perform a KS test on the two histograms
!			7. Print results of KS test to file, "ksOut.txt" and "dualKS.txt"
!			8. Set fluxCDF.nT to double its value. fluxCDF.nT gives the WE step that the next ks test is performed at.
!---------------------------------------------------------------------------------------------------------------------
*/
	FILE *KSFile, *dKSFile; //, *debugKS;
	int iLine, iBin, zeroCounts[3];
	double binStep, cdf1, cdf2, ksStat, dualCDF1, dualCDF2, dualKS;
	int nLines = nWE;
	double fluxVect[nWE];
	/* This section measures the KS by creating a histogram of the fluxes, using the histogram to measure the cdfs,
	then 
	*/
	//Clear previous binCounts
	for(iBin = 0; iBin < NFLUXBINS; iBin++){
		fluxCDF.binCounts[iBin] = 0;
		fluxCDF.oldCounts[iBin] = 0;
		fluxCDF.dualCounts[iBin] = 0;
		fluxCDF.oldDualCounts[iBin] = 0;
	}

	//load flux
	zeroCounts[0] = 0; //Necessary parameter for normalizing the non-zero CDF.
	zeroCounts[1] = 0;
	zeroCounts[2] = 0;
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
	if(fluxCDF.fluxMax > 0){
	for(iBin = 1; iBin <= NFLUXBINS; iBin++){
		fluxCDF.binDefs[iBin] = (double)iBin*binStep;
	}
	
	//Update Bin Counts
	for(iLine = (nLines/3) - 1; iLine < (2*nLines/3); iLine++){
		for(iBin=0;iBin<NFLUXBINS;iBin++){
			if(fluxVect[iLine]>=fluxCDF.binDefs[iBin] && fluxVect[iLine]<fluxCDF.binDefs[iBin+1]){
				fluxCDF.oldCounts[iBin]++;
				fluxCDF.oldDualCounts[iBin]++;
			}
		}
		if(fluxVect[iLine]==0){
					zeroCounts[0]++;
		}
		if(zeroCounts[0] > iLine){
			printf("statement");
		}
	}
	
	for(iLine = (2*nLines/3); iLine < nLines; iLine++){
		for(iBin=0;iBin<NFLUXBINS;iBin++){
			if(fluxVect[iLine]>=fluxCDF.binDefs[iBin] && fluxVect[iLine]<fluxCDF.binDefs[iBin+1]){
				fluxCDF.binCounts[iBin]++;
				fluxCDF.dualCounts[iBin]++;
			}
		}
				
		if(fluxVect[iLine]==0){
					zeroCounts[1]++;
		}		
	}
	
	//Determine which of the inspected CDFs has the least amount of zeros, then adjust the bincounts to remove those zeros.
	if(zeroCounts[0] >= zeroCounts[1]){
		zeroCounts[2] = zeroCounts[1];
	}
	if(zeroCounts[0]<zeroCounts[1]){
		zeroCounts[2] = zeroCounts[0];
	}
	fluxCDF.dualCounts[0] -= zeroCounts[2];
	fluxCDF.oldDualCounts[0] -= zeroCounts[2];
	//Perform KS Test
	cdf1 = 0;
	cdf2 = 0;
	dualCDF1 = 0;
	dualCDF2 = 0;
	ksStat = 0;
	dualKS = 0;
	for(iBin = 0; iBin < NFLUXBINS; iBin++){
		cdf1 += (double)fluxCDF.oldCounts[iBin]/(fluxCDF.nT/3);
		cdf2 += (double)fluxCDF.binCounts[iBin]/(fluxCDF.nT/3);
		dualCDF1 += (double)fluxCDF.oldDualCounts[iBin]/(fluxCDF.nT/3 - zeroCounts[2]); //Total counts used = (1/3 of time being measured) - number of zeros subtracted
		dualCDF2 += (double)fluxCDF.dualCounts[iBin]/(fluxCDF.nT/3 - zeroCounts[2]); 
		if(fabs(cdf1-cdf2)>ksStat){
			ksStat = fabs(cdf1-cdf2);
		}
		if(fabs(dualCDF1-dualCDF2)>dualKS){
			dualKS = fabs(dualCDF1-dualCDF2);
		}
	}
		if(zeroCounts[2] >= 1*nLines/3){
			ksStat = 2;
			dualKS = 2;
		}
		if(fabs(zeroCounts[1]-nLines/3) < NNZMIN || fabs(zeroCounts[0]-nLines/3) < NNZMIN){
			dualKS = 2;
		}
	}
	else{
		dualKS = 2;
		ksStat = 2;
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
	
	dKSFile = fopen("dualKS.txt","a");
	fprintf(dKSFile, "maxFlux = %E \n", fluxCDF.fluxMax);
	for(iBin = 0; iBin < NFLUXBINS; iBin++){
		fprintf(dKSFile, "%i ", fluxCDF.dualCounts[iBin]);
	}
	fprintf(dKSFile, "\n");
	for(iBin = 0; iBin < NFLUXBINS; iBin++){
		fprintf(dKSFile, "%i ", fluxCDF.oldDualCounts[iBin]);
	}
	fprintf(dKSFile, "\n dksStat: %E \n nT %i \n \n Zeros counted (removed): %i %i (%i)\n",dualKS, fluxCDF.nT, zeroCounts[0],zeroCounts[1],zeroCounts[2]);
	fclose(dKSFile);
	

	
	fluxCDF.ksStat = ksStat;
	fluxCDF.dualKS = dualKS;
	fluxCDF.nT *= 2;
	return;
}

void bin1Entropy(int nWE){
/*
!---------------------------------------------------------------------------------------------------------------------
!		Description:
!			Measures entropy and number of splits occurring in bin1 (bin next to flux bin at time of writing)
!		Notes:
!			As of 9/21/21 this function is not used for any data published by R. Taylor
!		Method:
!			Gathers number of splits from ending replicas - number of current replicas. Calculates entropy through
!			S = - sum(pi * log(pi)), except to keep loops minimized the calculation is rewritten as
!			S = - 1/wT * sum(wi * log(wi)) + log(wT), where pi = wi/wT, wT = sum(wi), and wi is the weight of a rep.
!			After calculation appends the timestep, number of splits, and entropy to a text file.
!			nSplits < 0 means that there will be merging rather than splitting in the next split merge step
!		
!---------------------------------------------------------------------------------------------------------------------
*/
	FILE *binFile;
	int nReps,nSplits,iReps;
	double entropy, relativeEntropy, relativeWeightTotal, iWeight;
	binFile = fopen("bin1.txt","a");
	
	nReps = Reps.binContentsMax[1];
	relativeWeightTotal = 0;
	relativeEntropy = 0;
	for(iReps = 0; iReps < nReps; iReps++){
		iWeight = Reps.weights[Reps.binContents[iReps][1]];
		relativeWeightTotal+= iWeight;
		relativeEntropy+= - iWeight * log(iWeight);
	}
	entropy = relativeEntropy/relativeWeightTotal + log(relativeWeightTotal);
	nSplits = paramsWe.repsPerBin - nReps;
	fprintf(binFile,"%i, %i, %E\n",nWE, nSplits,entropy);
	fclose(binFile);
}

void evacuatedParticles(){
/*
!---------------------------------------------------------------------------------------------------------------------
!		Description:
!			Measures particle locations + replica weights for systems in flux bin
!		Notes:
!			As of 9/21/21 this function is not used for any data published by R. Taylor.
!		Method:
!			For each sim in the flux bin, write its weight + monomer locations to 1 file and weight + dimer locations
!			to another file.
!				1. Identify sims in flux bin
!				2. For each sim, write its weight to monomer + dimer file
!				3. For both lists of molecules, loop through all live molecules, writing x+y pos to same line of files
!				4. After writing positions to file, write new line character
!---------------------------------------------------------------------------------------------------------------------
*/
	FILE *monomersFile, *dimersFile;
	int iBC, nList, iSim, nMol;
	int nFlux = Reps.binContentsMax[paramsWe.fluxBin];
	if(nFlux>0){
	monomersFile = fopen("monomerLocs.txt","a");
	dimersFile = fopen("dimerLocs.txt","a");
		for(iBC = nFlux -1; iBC >=0; iBC--){ //Loop through the replicas
			iSim = Reps.binContents[iBC][paramsWe.fluxBin];
			fprintf(monomersFile,"%e, ", Reps.weights[iSim]);
			fprintf(dimersFile,"%e, ",Reps.weights[iSim]);
				for(nList = 0; nList < Reps.sims[iSim]->mols->nlist;nList++){
					if(nList == 0){
						for(nMol = 0; nMol < Reps.sims[iSim]->mols->nl[nList];nMol++){
							fprintf(monomersFile,"%e, %e, ",Reps.sims[iSim]->mols->live[nList][nMol]->pos[0],Reps.sims[iSim]->mols->live[nList][nMol]->pos[1]);
						}
						fprintf(monomersFile,"\n");
					}
					if(nList ==1){
						for(nMol = 0; nMol < Reps.sims[iSim]->mols->nl[nList];nMol++){
							fprintf(dimersFile,"%e, %e, ",Reps.sims[iSim]->mols->live[nList][nMol]->pos[0],Reps.sims[iSim]->mols->live[nList][nMol]->pos[1]);
						}
						fprintf(dimersFile,"\n");
					}
				}
		}
		fclose(monomersFile);
		fclose(dimersFile);
	}
	
	
}

void saveWE(){
	
	/*
!---------------------------------------------------------------------------------------------------------------------
!		Description:
!			Measures particle locations + replica weights to save WE state
!		Method:
!			For each sim write its weight + monomer locations to a file to be accessed later
!				1. Loop through all active sims
!				2. Record replica weight
!				3. Write monomer X + Y locations
!				4. Write dimer X + Y locations
!---------------------------------------------------------------------------------------------------------------------
*/
	FILE *stateFile;
	int iSim, nList, nMol;
	
	stateFile = fopen("savestate.txt","w");
	fprintf(stateFile,"iSimMax %i \n",Reps.iSimMax);
	for(iSim =0; iSim < Reps.iSimMax; iSim++){
		fprintf(stateFile,"Weight %e \n",Reps.weights[iSim]);
		for(nList = 0; nList < Reps.sims[iSim]->mols->nlist;nList++){
			if(nList == 0){
				fprintf(stateFile,"Monomers %i\n",Reps.sims[iSim]->mols->nl[nList]);
				for(nMol = 0; nMol < Reps.sims[iSim]->mols->nl[nList];nMol++){
						fprintf(stateFile,"%le, %le \n",Reps.sims[iSim]->mols->live[nList][nMol]->pos[0],Reps.sims[iSim]->mols->live[nList][nMol]->pos[1]);
					}
				}
			if(nList ==1){
				fprintf(stateFile,"Dimers %i\n",Reps.sims[iSim]->mols->nl[nList]);
				for(nMol = 0; nMol < Reps.sims[iSim]->mols->nl[nList];nMol++){
						fprintf(stateFile,"%le, %le \n",Reps.sims[iSim]->mols->live[nList][nMol]->pos[0],Reps.sims[iSim]->mols->live[nList][nMol]->pos[1]);
					}
				}
			}
	}
	fclose(stateFile);
}

void loadWE(){
/*
!---------------------------------------------------------------------------------------------------------------------
!		Description:
!			Uses previous savestate file to initialize the WE simulation. 
!		Notes:
!			This only executes if a command-line argument, the load bit, is set to 1, though the function itself does
!			not reference it.
!		Method:
!			1. Reads the savestate file for each replica. Creates a new empty replica with no molecules inside
!			2. Gives the empty sim the WE weight given by the savestate
!			3. Adds molecules to the sim at the locations given by the savestate.
!---------------------------------------------------------------------------------------------------------------------
*/
	FILE *stateFile;
	char tmpStr[32];
	int iSim, nMonomers, nDimers,iMol;
	double molPos[2];
	char const* monomer = "A";
	char const* dimer = "B";
	char **stateList;
	stateList = (char**) calloc(2,sizeof(char*));
	stateList[0] = (char*) calloc(256, sizeof(char));
	stateList[1] = (char*) calloc(256, sizeof(char));
	strcpy(stateList[0],monomer);
	strcpy(stateList[1],dimer);
	
	stateFile = fopen("savestate.txt","r");
	
	fscanf(stateFile,"%s %i",tmpStr, &Reps.iSimMax);
	for(iSim = 0; iSim < Reps.iSimMax;iSim++){
		buildEmptySim(iSim);
		fscanf(stateFile,"%s %le",tmpStr, &Reps.weights[iSim]);
		fscanf(stateFile,"%s %i",tmpStr, &nMonomers);
		for(iMol = 0; iMol < nMonomers; iMol++){
			fscanf(stateFile,"%le, %le",&molPos[0],&molPos[1]);
			smolAddSolutionMolecules(Reps.sims[iSim],monomer, 1,molPos,molPos);
			smolUpdateSim(Reps.sims[iSim]);
		}
		fscanf(stateFile,"%s %i",tmpStr,&nDimers);
		for(iMol = 0; iMol < nDimers; iMol++){
			fscanf(stateFile,"%le, %le",&molPos[0],&molPos[1]);
			smolAddSolutionMolecules(Reps.sims[iSim],dimer,1,molPos,molPos);
		}
		Reps.binLocs[iSim]=findBin(Reps.sims[iSim]);
	}
}

int main(int argc, char *argv[]){
	/*
!---------------------------------------------------------------------------------------------------------------------
!		Description:
!			Main WE function. Intakes command-line arguments, reads the parameters files, initializes the sims, 
!			writes to files. Additionally tracks things like monomerization fraction, time elapsed.
!		Notes:
!			argv1: ending simfile (Somewhat of a misnomer, gives weight distribution across bins at given timesteps)
!			argv2: flux file. Appended to every WE timestep with the flux measured in that timestep
!			argv3: Seed/error file. Outputs the ISEED used + miscellaneous errors that I wanted to output here.
!			argv4: Save / replace RNG bit. If this is 1 then the RNG seed is saved / ISEED does not increment. If 0
!					then it increments according to twister.c
!			argv5: Execution time file. Outputs, in order: 1. Time spent creating initial distribution 2. All time
!					spent in the WE splitting/merging 3. All time spent in Smoldyn dynamics 4. Total time of 
!					WE/Dyanamics 5. Time spent creating the savestate
!			argv6: Load sim bit. If 1, then savestate.txt is referenced to load the previous savestate, if 0 then
!					this initializes a new WE sim.
!		Method:
!			1. Initialize everything
!			2. begin WE simulation loop:
!				a. Records fluxes
!				b. Records KS statistics
!				c. Execute WE splitting and merging
!				d. Execute Smoldyn dynamics for each replica. Update bin locations. Update continuing monomerization
!					fractions
!			3. Record all to files.
!---------------------------------------------------------------------------------------------------------------------
*/
	//argv 1: ending simfile, argv2: flux file, argv3: seed / error file, argv4: save / replace rng bit argv5: Execution time file argv6: Load Sim Bit (1 = load savestate.txt)
	//dynamics params: dt, L, R, D, N
	//WE Params: tau, mTarg, tauMax, nBins, ((flux bin))
	int tauMax, rngBit, iBin, nWE, iSim, iBCM, nanCheck, firstNAN, iDimer,iBinContents,iClockInit,mCounts[3],tauInit,loadBit,nWEstart; //tauQuarter omitted
	double fluxAtStep, binWeight,mCountsWeighted[3], clockDouble[5];
	clock_t start[5], stop[5]; //initialDistTime, splitMergeTime, dynamicsTime, totalTime, saveTime, this also corresponds to the order written in the output file
	char fileReader;
	//Load simulation / WE parameters from outside files
	FILE *DEFile, *WEFile, *CorralsFile, *FLFile, *SIMFile, *errFile, *clockFile, *mCountsFile, *structStoreFile, *structNANFile, *debugFile, *mCountsWeightedFile, *binTimeSeriesFile;

	char *fluxFileStr;
	if(1){ //the only reason this exists is so that i can minimize this part while editing
	fluxFileStr = argv[2];
	DEFile = fopen("dynamicsParams.txt","r");
	WEFile = fopen("WEParams.txt","r");
	CorralsFile = fopen("corralsParams.txt","r");
	errFile = fopen(argv[3], "w");
	
	getParams(DEFile, WEFile,CorralsFile);
	
	loadBinDefs();
		
	//Defining variables for continuous dimer counting, gives more data points for comparing monomerization fraction
	
	//These files will each give the total monomers, dimers, and either all "time" elapsed for all sims or all weight measured across
	for(iDimer = 0; iDimer < 3; iDimer++){
		mCounts[iDimer] = 0;
		mCountsWeighted[iDimer] = 0;
	}
	for(iClockInit = 0; iClockInit<=4; iClockInit++){
		clockDouble[iClockInit]=0;
	}
	tauMax = paramsWe.tauMax;
		tauInit = paramsWe.tauMax;
	if(DEBUGGING){
		debugFile = fopen("Debug.txt","a");
		fprintf(debugFile,"Tau loops + flux vector made \n");
		fclose(debugFile);
	}
	//Set + record new RNG seed, set up initial distribution
	rngBit = atoi(argv[4]);
	loadBit = atoi(argv[6]);
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
		if(~loadBit){
			initialDist(paramsWe.nInit);
			nWEstart = 0;
		}
		if(loadBit){
			loadWE();
			nWEstart = 1;
			FLFile = fopen(fluxFileStr,"r");
			for(fileReader = getc(FLFile);fileReader!=EOF;fileReader=getc(FLFile)){
				if(fileReader =='\n'){
					nWEstart +=1;
				}
			}
			if(nWEstart>=tauMax){
				tauMax *=2;
			}
		}
	stop[0] = clock();
	clockDouble[0]+=(double)(stop[0]-start[0])/CLOCKS_PER_SEC;
	if(DEBUGGING){
	debugFile = fopen("Debug.txt","a");
	fprintf(debugFile, "Initial Distribution Made \n");
	fclose(debugFile);
	}
	//Set other key initial values in replicas struct
	for(iBin = 0; iBin < Reps.nBins; iBin++){
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
	}
	
	/*As of right now, I want to set up some structure to measure how monomerization fraction 
	and distribution of bin weights changes over time. For each bin, I am going to create an 
	output file that records the monomerization fraction, weight inside the bins, at each tau step
	That is many files, and the amount of files changes based on the number of bins*/
	
	//Simulation Loop
	for(nWE = nWEstart-1; nWE < tauMax; nWE++){
		if(nWE > nWEstart || !loadBit){
		//Print current tau step to stdout for interactive debugging
		if(DEBUGGING){
		debugFile = fopen("Debug.txt","a");
		fprintf(debugFile,"Tau Step: %i \n", nWE);
		fclose(debugFile);
		}
		//Flux Recording
		if(paramsWe.fluxBin >= 0){
			evacuatedParticles();
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
			
			//Makeshift bug fix when paramsWe.replaceFluxSims==1
			if(paramsWe.replaceFluxSims==1){
						for( iBin = 0; iBin < Reps.nBins; iBin++){
			Reps.binContentsMax[iBin] = 0;
		}
		for(iSim = 0; iSim < Reps.iSimMax; iSim++){

			Reps.binLocs[iSim] = findBin(Reps.sims[iSim]);
			Reps.binContents[Reps.binContentsMax[Reps.binLocs[iSim]]][Reps.binLocs[iSim]] = iSim;
			Reps.binContentsMax[Reps.binLocs[iSim]]++;
	}
			}
		
		//Single Step Time Measurement
		
			start[1] = clock();
		
		
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
			if(!FIXEDTIME){
			if(fluxCDF.dualKS >= DKSCRITICAL || fluxCDF.ksStat >= KSCRITICAL){
				tauMax += tauMax;
			}
			
			if(fluxCDF.dualKS < DKSCRITICAL && fluxCDF.ksStat < KSCRITICAL){
				nWE = tauMax;
			}
			}
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
			stop[1] = clock();
		}
		//Saving
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
		
		clockDouble[1]+=(double)(stop[1]-start[1])/CLOCKS_PER_SEC;
		start[4] = clock();
		if(nWE % 10 == 0){	
		saveWE();
		}
		stop[4] = clock();
		clockDouble[4] += (double)(stop[4]-start[4])/CLOCKS_PER_SEC;
		if(clockDouble[1]+clockDouble[4]+clockDouble[2]>250000){
			nWE = tauMax;
		}
		if(DEBUGGING){
			debugFile = fopen("Debug.txt","a");
			fprintf(debugFile,"NANs checked and recorded. iSimMax = %i \n Dynamics Starting \n", Reps.iSimMax);
			fclose(debugFile);
		}
		
		//Dynamics and rebuilding BCM
		for( iBin = 0; iBin < Reps.nBins; iBin++){
			Reps.binContentsMax[iBin] = 0;
		}
		for(iSim = 0; iSim < Reps.iSimMax; iSim++){
			/*if(DEBUGGING){
				debugFile = fopen("Debug.txt","a");
				fprintf(debugFile,"Sim %i dynamics starting  \n", iSim);
				fclose(debugFile);
						 }*/
			start[2] = clock();
			dynamicsEngine(Reps.sims[iSim]);
			stop[2] = clock();
			clockDouble[2]+=(double)(stop[2]-start[2])/CLOCKS_PER_SEC;
			Reps.binLocs[iSim] = findBin(Reps.sims[iSim]);
			Reps.binContents[Reps.binContentsMax[Reps.binLocs[iSim]]][Reps.binLocs[iSim]] = iSim;
			Reps.binContentsMax[Reps.binLocs[iSim]]++;
			//Dimerization Fraction Recording
			if(!MONOFRACEACHDT){
 			if(nWE > tauInit/2){
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
		binTimeSeriesFile = fopen("timeSeries.txt","a");
		double mCountTemp, weightTemp,simTimeTemp;
		for(iBin = 0; iBin < Reps.nBins; iBin++){
			weightTemp = 0;
			mCountTemp = 0;
			simTimeTemp = nWE;
			for(iSim=0; iSim < Reps.binContentsMax[iBin];iSim++){
				weightTemp += Reps.weights[Reps.binContents[iSim][iBin]];
				mCountTemp += (double) Reps.sims[Reps.binContents[iSim][iBin]]->mols->nl[0] * Reps.weights[Reps.binContents[iSim][iBin]];
			}
			fprintf(binTimeSeriesFile,"%d %f %f %d \n", iBin, weightTemp, mCountTemp/(weightTemp*paramsDe.nPart), nWE);
		}
		
		fclose(binTimeSeriesFile);
		
		
		
		
		
		
		
		
		
		
		
		
		if(tauMax >= 1000){
		if( (nWE+1) % (tauMax/1000) == 0 || nWE < 1000){
			bin1Entropy(nWE);
		SIMFile = fopen(argv[1],"a");	
		fprintf(SIMFile,"tau = %i\n",nWE);
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
			bin1Entropy(nWE);
		SIMFile = fopen(argv[1],"a");
		fprintf(SIMFile,"tau = %i\n",nWE);
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
		}
		
		if(DEBUGGING){
			debugFile = fopen("Debug.txt","a");
			fprintf(debugFile,"Dynamics Finished  \n Recording mCounts \n");
			fclose(debugFile);
		}
		
		if(nWE > tauInit/2){
			//mCountsFile = fopen("mCounts.txt","a");
			//fprintf(mCountsFile, "%i, %i, %i \n", mCounts[0], mCounts[1], mCounts[2]);
			//fclose(mCountsFile);
			mCountsWeightedFile = fopen("mCountsWeighted.txt","a");
			fprintf(mCountsWeightedFile, "%f, %f, %f \n", mCountsWeighted[0], mCountsWeighted[1],mCountsWeighted[2]);
			fclose(mCountsWeightedFile);
		}
	}
	
	stop[3] = clock();
	clockDouble[3]+=(double)(stop[3]-start[3])/CLOCKS_PER_SEC;
	
	//Time Recording
	clockFile = fopen(argv[5],"a");
	fprintf(clockFile,"%E \n%E \n%E \n%E \n%E \n",clockDouble[0],clockDouble[1],clockDouble[2],clockDouble[3],clockDouble[4]);
	fclose(clockFile);
	
	//Free Memory and finish
	freeAllSims();
	return 0;
}