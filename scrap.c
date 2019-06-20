/*Functions to build: 
findOrderParameter
orderParameterToBin
LoadBinDefs
dimerEquilDist.

Need to modify globals to account for:
Order parameters as they account for bins

The goal for setting up the order parameters / bins should be to tile the order parameter space
Specify min / max order parameters:

Need to clarify number of order parameters somewhere in the code. say "n" order parameters

For each bin, we should specify a n-d box that the bin covers
e.g. [0,1][0,1] for a box of area 1 with corner at the origin (2 order parameters).

For now: assume that there will only ever be two order parameters: dimer count and number inside the roi.

Should have some sort of file that tells us important bin information. I'll make a "binParams.txt" file.
*/
struct BinDefinitions{
	double binDefArray[4*NBINSMAX];
	int nBins;
	int currentDims;
	int customBins;
}

struct BinDefinitions binDefs;

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
	int useDimerDist;
} ;


void loadBinDefs(){
	FILE *binDefinitions, *binParams;
	char tmpChar[32];
	int iBin, iNbin;
	binParams = fopen("binParams.txt","r");
	binDefinitions = fopen("binDefinitions.txt","r");
	fscanf(binParams,"%s %i", tmpChar, &binDefs.customBins);
	fscanf(binParams,"%s %i", tmpChar, &binDefs.currentDims);
	fscanf(binParams,"%s %i",tmpChar, &binDefs.nBins);
	fclose(binParams);
	if(binDefs.customBins == 1){
	for(int iBinDef = 0; iBinDef < (binDefs.currentDims*binDefs.currentDims*binDefs.nBins); iBinDef++){
		fscanf(binDefinitions,"%i", &binDefs.binDefArray[iBinDef]);
	}
	}
}

int findOrderParameter(simptr currentSim, int orderParamIndex){
	int nInBin, nMol, nMolList;
	double molX, molY;
	orderParam = 0;
	
	if(orderParamIndex == ROIORDERPARAM){
		for(nMolList = 0; nMolList < currentSim->mols->nlist; nMolList++){ 
		for(nMol = 0; nMol < currentSim->mols->nl[nMolList]; nMol++){
		molX = currentSim->mols->live[nMolList][nMol]->pos[0]; //first index changes based on how we org mol lists
		molY = currentSim->mols->live[nMolList][nMol]->pos[1];
		if((molX*molX+ molY*molY) < paramsDe.roiR*paramsDe.roiR){
			nInBin++;
		}
		}
		}
	}
	
	if(orderParamIndex == DIMERORDERPARAM){
		orderParam = currentSim->mols->nl[1];
	}
	
	return orderParam;
}

int findBin(simptr currentSim){
	int ordParam1, ordParam2, iBin, binOut;
	
	ordParam1 = findOrderParameter(currentSim, ROIORDERPARAM);
	ordParam2 = findOrderParameter(currentSim, DIMERORDERPARAM);
	
	if(binDefs.customBins == 1){
	for(iBin = 0; iBin < nBins; iBin++){
		if (binDefs.currentDims==2){
		if((ordParam1 >= binDefs[4*iBin]) && (ordParam1 <= binDefs[4*iBin+1]) && (ordParam2 >= binDefs[4*iBin+2]) && (ordParam2 <= binDefs[4*iBin+3])){
			binOut = iBin;
			iBin = nBins;
		}
		}
		if(binDefs.currentDims==1){
			if((ordParam1 >= binDefs[2*iBin])&&(ordParam1>=binDefs[]))
		}
	}
	}
	else{
		if(ROBINS){
			binOut = ordParam1;
		}
		if(!ROBINS){
			binOut = ordParam1;
		}
	}
	
	return binOut;
}

void initialDist(int nInit){
	int jSim, iDimer, nDimers,tmpDimer;
	double lowBounds[] = {-paramsDe.worldLength/2,-paramsDe.worldLength/2};
	double highBounds[] = {paramsDe.worldLength/2, paramsDe.worldLength/2};
	double botLeftCornerRect[] = {-paramsDe.worldLength/2, -paramsDe.worldLength/2, paramsDe.worldLength};
	double topRightCornerRect[] = {paramsDe.worldLength/2, paramsDe.worldLength/2, -paramsDe.worldLength};
	double roiParams[] = {0.0, 0.0, paramsDe.roiR, 30};
	double insideRoi[] = {0.0, 0.0};
	char const* monomer = "A";
	char const* dimer = "B";
	char const** dimerptr = &dimer;
	enum MolecState outputStates[2];
	char **unbindProducts;
	unbindProducts = (char**) calloc(2,sizeof(char*));
	unbindProducts[0] = (char*) calloc(256,sizeof(char));
	unbindProducts[1] = (char*) calloc(256,sizeof(char));
	strcpy(unbindProducts[0],monomer);
	strcpy(unbindProducts[1],monomer);
	outputStates[0] = MSsoln;
	outputStates[1] = MSsoln;
	
	FILE *dimerCDFfile;
	double dimerDist[paramsDe.nPart/2], tmpRand;
	tmpRand = RAND;
	nDimers = 0;
	if(paramsDE.useDimerDist == 1){
		dimerCDFfile = fopen("dimerDistrbution.txt")
			for(iDimer = 0; iDimer < paramsDe.nPart/2;iDimer++){
				fscanf(dimerCDFfile,"%i %lf",&tmpDimer, dimerDist[iDimer]);
				if(tmpRand < dimerDist[iDimer]){
					nDimers = tmpDimer;
				}
			}
		fclose(dimerDistFile);
		
	}
	
	for(jSim = 0; jSim < nInit; jSim++){
		//Molecules + BCs
		Reps.sims[jSim] = smolNewSim(2, lowBounds,highBounds);
		smolSetRandomSeed(Reps.sims[jSim],genrand_int31());
		smolSetGraphicsParams(Reps.sims[jSim], "none", 1, 0);
		smolSetSimTimes(Reps.sims[jSim],0,10000,paramsDe.dt);
		
		smolAddMolList(Reps.sims[jSim],"mList");
		smolAddSpecies(Reps.sims[jSim],monomer,"mList");
		smolSetSpeciesMobility(Reps.sims[jSim],monomer,MSall, paramsDe.difM, NULL, NULL);
		smolSetMaxMolecules(Reps.sims[jSim],2*paramsDe.nPart);
		
		if(paramsDe.reactBit > 0){
			smolAddMolList(Reps.sims[jSim],"dList");
			smolAddSpecies(Reps.sims[jSim], dimer, "dList");
			smolSetSpeciesMobility(Reps.sims[jSim],dimer, MSall, paramsDe.difD, NULL, NULL);
			smolAddReaction(Reps.sims[jSim], "binding", monomer, MSsoln, monomer, MSsoln, 1, dimerptr, &outputStates[0], -1);
			smolSetReactionRate(Reps.sims[jSim], "binding", paramsDe.bindR, 1); //This line allows us to set the rate by the binding radius rather than kOn
			smolAddReaction(Reps.sims[jSim], "unbinding", dimer, MSsoln, NULL , MSnone, 2,(const char**) unbindProducts, outputStates, paramsDe.unbindK);
		}
		
		smolAddSolutionMolecules(Reps.sims[jSim], monomer, paramsDe.nPart-2*nDimers, lowBounds, highBounds);
		if(paramsDe.reactBit>0){
			smolAddSolutionMolecules(Reps.sims[jSim],dimer,nDimers,lowBounds,highBounds);
		}
		smolAddSurface(Reps.sims[jSim], "bounds");
		smolAddPanel(Reps.sims[jSim], "bounds", PSrect, NULL, "-x", topRightCornerRect);
		smolAddPanel(Reps.sims[jSim], "bounds", PSrect, NULL, "-y", botLeftCornerRect);
		smolAddPanel(Reps.sims[jSim], "bounds", PSrect, NULL, "+x", botLeftCornerRect);
		smolAddPanel(Reps.sims[jSim], "bounds", PSrect, NULL, "+y", topRightCornerRect);
		smolSetSurfaceAction(Reps.sims[jSim], "bounds", PFboth, "all", MSall, SAreflect);
		
		if(STOPCOMMAND){
		smolAddCommandFromString(Reps.sims[jSim], "e ifincmpt all = 0 roiComp stop");
		}
		
		//ROI Surface + compartment
		smolAddSurface(Reps.sims[jSim], "roi");
		smolAddPanel(Reps.sims[jSim],"roi", PSsph, NULL, "+0", roiParams);
		smolSetSurfaceAction(Reps.sims[jSim], "roi", PFboth, "all", MSall, SAtrans);
		smolAddCompartment(Reps.sims[jSim],"roiComp");
		smolAddCompartmentSurface(Reps.sims[jSim],"roiComp","roi");
		smolAddCompartmentPoint(Reps.sims[jSim],"roiComp",insideRoi);
		smolUpdateSim(Reps.sims[jSim]);
		
		//Reps.sims[jSim]->logfile = nulldev;
	}
}