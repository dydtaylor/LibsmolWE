void initialDist(int nInit){
	int jSim;
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
		
		smolAddSolutionMolecules(Reps.sims[jSim], monomer, paramsDe.nPart, lowBounds, highBounds);
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


void dynamicsEngine(simptr currentSim, int mCounts[], double mCountsWeighted[], int nWE, int iSim){
	smolUpdateSim(currentSim);
	for(int i = 0; i < paramsWe.tau; i++){
	smolRunTimeStep(currentSim);
		if(MONOFRACEACHDT){
			if(nWE > paramsWe.tauMax /2){
				mCounts[0] += smolGetMoleculeCount(currentSim, "A", MSall);
				mCounts[1] += smolGetMoleculeCount(currentSim, "B", MSall);
				mCounts[2] += 1;
				mCountsWeighted[0] += (double) smolGetMoleculeCount(currentSim,"A", MSall) * Reps.weights[iSim];
				mCountsWeighted[1] += (double) smolGetMoleculeCount(currentSim,"B", MSall) * Reps.weights[iSim];
				mCountsWeighted[2] += Reps.weights[iSim];
	}
	}
	}
}

int findOrderParameter(simptr currentSim, int orderParamIndex){
	int nMol, nMolList, orderParam;
	double molX, molY;
	orderParam = 0;
	
	if(orderParamIndex == ROIORDERPARAM){
		for(nMolList = 0; nMolList < currentSim->mols->nlist; nMolList++){ 
		for(nMol = 0; nMol < currentSim->mols->nl[nMolList]; nMol++){
		molX = currentSim->mols->live[nMolList][nMol]->pos[0]; //first index changes based on how we org mol lists
		molY = currentSim->mols->live[nMolList][nMol]->pos[1];
		if((molX*molX+ molY*molY) < paramsDe.roiR*paramsDe.roiR){
			orderParam++;
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
	for(iBin = 0; iBin < binDefs.nBins; iBin++){
		if (binDefs.currentDims==2){
		if((ordParam1 >= binDefs.binDefArray[4*iBin]) && (ordParam1 < binDefs.binDefArray[4*iBin+1]) && (ordParam2 >= binDefs.binDefArray[4*iBin+2]) && (ordParam2 < binDefs.binDefArray[4*iBin+3])){
			binOut = iBin;
			iBin = binDefs.nBins;
		}
		}
		if(binDefs.currentDims==1 && ROBINS){
			if((ordParam1 >= binDefs.binDefArray[2*iBin])&&(ordParam1<binDefs.binDefArray[2*iBin+1])){
				binOut=iBin;
				iBin = binDefs.nBins;
			}
		}
		if(binDefs.currentDims==1&& (!ROBINS)){
			if((ordParam2 >= binDefs.binDefArray[2*iBin])&&(ordParam2<binDefs.binDefArray[2*iBin+1])){
				binOut=iBin;
				iBin = binDefs.nBins;
			}
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

void freeAllSims(){
	int jSim;
	for(jSim = 0; jSim <=Reps.iSimMax; jSim++){
	smolFreeSim(Reps.sims[jSim]);
	}
}

void copySim1(int simIn, int simOut){

	/*
	Takes simIn, reads molecule locations from the sim, then makes a new sim with molecules at the same locations. Most parameters / lines of code are identical to the sim made through initialDist
	*/
	
	int nList, nMol;
	//FILE *nulldev = fopen(NULLDEVICE, "w");
	
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
	char **unbindProducts, **stateList;
	unbindProducts = (char**) calloc(2,sizeof(char*));
	stateList = (char**) calloc(2,sizeof(char*));
	unbindProducts[0] = (char*) calloc(256,sizeof(char));
	unbindProducts[1] = (char*) calloc(256,sizeof(char));
	stateList[0] = (char*) calloc(256, sizeof(char));
	stateList[1] = (char*) calloc(256, sizeof(char));
	strcpy(unbindProducts[0],monomer);
	strcpy(unbindProducts[1],monomer);
	strcpy(stateList[0],monomer);
	strcpy(stateList[1],dimer);
	outputStates[0] = MSsoln;
	outputStates[1] = MSsoln;
	
	smolFreeSim(Reps.sims[simOut]);
	Reps.sims[simOut] = smolNewSim(2, lowBounds, highBounds);
	smolSetRandomSeed(Reps.sims[simOut],genrand_int31());
	smolSetGraphicsParams(Reps.sims[simOut], "none", 1, 0);
	smolSetSimTimes(Reps.sims[simOut], 0, 10000, paramsDe.dt); //Reps.sims[simIn]->time is another option
	smolAddMolList(Reps.sims[simOut],"mList");
	smolAddSpecies(Reps.sims[simOut], "A", "mList");
	smolSetSpeciesMobility(Reps.sims[simOut], "A", MSall, paramsDe.difM, NULL, NULL);
	smolSetMaxMolecules(Reps.sims[simOut],2*paramsDe.nPart);
	
		if(paramsDe.reactBit > 0){
			smolAddMolList(Reps.sims[simOut],"dList");
			smolAddSpecies(Reps.sims[simOut], dimer, "dList");
			smolSetSpeciesMobility(Reps.sims[simOut],dimer, MSall, paramsDe.difD, NULL, NULL);
			smolAddReaction(Reps.sims[simOut], "binding", monomer, MSsoln, monomer, MSsoln, 1, dimerptr, &outputStates[0], -1);
			smolSetReactionRate(Reps.sims[simOut], "binding", paramsDe.bindR, 1); //This line allows us to set the rate by the binding radius rather than kOn
			smolAddReaction(Reps.sims[simOut], "unbinding", dimer, MSsoln, NULL , MSnone, 2,(const char**) unbindProducts, outputStates, paramsDe.unbindK);
	}
	
	smolAddSurface(Reps.sims[simOut],"bounds");
	smolAddPanel(Reps.sims[simOut], "bounds", PSrect, NULL, "-x", topRightCornerRect);
	smolAddPanel(Reps.sims[simOut], "bounds", PSrect, NULL, "-y", botLeftCornerRect);
	smolAddPanel(Reps.sims[simOut], "bounds", PSrect, NULL, "+x", botLeftCornerRect);
	smolAddPanel(Reps.sims[simOut], "bounds", PSrect, NULL, "+y", topRightCornerRect);
	smolSetSurfaceAction(Reps.sims[simOut], "bounds", PFboth, "all", MSall, SAreflect);
		
	smolAddSurface(Reps.sims[simOut], "roi");
	smolAddPanel(Reps.sims[simOut],"roi", PSsph, NULL, "+0", roiParams);
	smolSetSurfaceAction(Reps.sims[simOut], "roi", PFboth, "all", MSall, SAtrans);
	smolAddCompartment(Reps.sims[simOut],"roiComp");
	smolAddCompartmentSurface(Reps.sims[simOut],"roiComp","roi");
	smolAddCompartmentPoint(Reps.sims[simOut],"roiComp",insideRoi);
	

	if(STOPCOMMAND){
	smolAddCommandFromString(Reps.sims[simOut], "e ifincmpt all = 0 roiComp stop");
	}
	//The following for loop is the major difference between this and the initial configuration simulation building.
	
	
	for(nList = 0; nList < Reps.sims[simIn]->mols->nlist; nList++){
		for(nMol = 0; nMol < Reps.sims[simIn]->mols->nl[nList];nMol++){
			smolAddSolutionMolecules(Reps.sims[simOut], stateList[nList], 1, Reps.sims[simIn]->mols->live[nList][nMol]->pos,Reps.sims[simIn]->mols->live[nList][nMol]->pos);
	}
	}
	smolUpdateSim(Reps.sims[simOut]);
	smolUpdateSim(Reps.sims[simIn]);
	//Print statement to check functionality. Note molecule indices get changed but everything else preserved
	/*
	for(nList = 0; nList < Reps.sims[simIn]->mols->nlist; nList++){
		for(nMol = 0; nMol < Reps.sims[simIn]->mols->nl[nList];nMol++){
		printf("Molecule %i Coords: %f %f \n",nMol, Reps.sims[simOut]->mols->live[nList][nMol]->pos[0],Reps.sims[simOut]->mols->live[nList][nMol]->pos[1]);
		printf("Molecule %i Coords: %f %f \n",nMol, Reps.sims[simIn]->mols->live[nList][nMol]->pos[0],Reps.sims[simIn]->mols->live[nList][nMol]->pos[1]);
		}
	}
	*/
}

