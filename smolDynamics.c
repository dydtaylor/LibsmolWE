/*void simSetup(int simOut){

	Reps.sims[simOut] = smolNewSim(2, lowBounds, highBounds);
	smolSetRandomSeed(Reps.sims[simOut],genrand_int31());
	smolSetGraphicsParams(Reps.sims[simOut], "none", 1, 0);
	smolSetSimTimes(Reps.sims[simOut], 0, SMOLTIMEMAX,paramsDe.dt);
	smolAddSpecies(Reps.sims[simOut],"A",NULL);
	smolSetSpeciesMobility(Reps.sims[simOut],"A",MSall, paramsDe.difC, 0, 0);
	smolSetMaxMolecules(Reps.sims[simOut],paramsDe.nPart);
	smolAddSurface(Reps.sims[simOut], "bounds");
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
	smolAddCommandFromString(Reps.sims[simOut], "e ifincmpt A = 0 roiComp stop");
	smolUpdateSim(Reps.sims[simOut]);
}

void initialDist(int nInit){
	int jSim;
	for(jSim = 0; jSim < nInit; jSim++){
		//Molecules + BCs
		simSetup(jSim);
		smolAddSolutionMolecules(Reps.sims[jSim], "A", paramsDe.nPart, NULL, NULL);
		smolUpdateSim(Reps.sims[jSim]);
		//Reps.sims[jSim]->logfile = nulldev;
	}
}*/
void initialDist(int nInit){
	int jSim;
	double lowBounds[] = {-paramsDe.worldLength/2,-paramsDe.worldLength/2};
	double highBounds[] = {paramsDe.worldLength/2, paramsDe.worldLength/2};
	double botLeftCornerRect[] = {-paramsDe.worldLength/2, -paramsDe.worldLength/2, paramsDe.worldLength};
	double topRightCornerRect[] = {paramsDe.worldLength/2, paramsDe.worldLength/2, -paramsDe.worldLength};
	double roiParams[] = {0.0, 0.0, 1, 30};
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
		smolSetGraphicsParams(Reps.sims[jSim], "opengl", 1, 0);
		smolSetSimTimes(Reps.sims[jSim],0,10000,paramsDe.dt);
		
		smolAddMolList(Reps.sims[jSim],"mList");
		smolAddSpecies(Reps.sims[jSim],monomer,"mList");
		smolSetSpeciesMobility(Reps.sims[jSim],monomer,MSall, paramsDe.difM, NULL, NULL);
		smolSetMaxMolecules(Reps.sims[jSim],paramsDe.nPart);
		smolAddSolutionMolecules(Reps.sims[jSim], monomer, paramsDe.nPart, lowBounds, highBounds);
		smolAddSurface(Reps.sims[jSim], "bounds");
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
		
		if(paramsDe.reactBit > 0){
			smolAddMolList(Reps.sims[jSim],"dList");
			smolAddSpecies(Reps.sims[jSim], dimer, "dList");
			smolSetSpeciesMobility(Reps.sims[jSim],dimer, MSall, paramsDe.difD, NULL, NULL);
			smolAddReaction(Reps.sims[jSim], "binding", monomer, MSsoln, monomer, MSsoln, 1, dimerptr, &outputStates[0], -1);
			smolSetReactionRate(Reps.sims[jSim], "binding", paramsDe.bindR, 1); //This line allows us to set the rate by the binding radius rather than kOn
			smolAddReaction(Reps.sims[jSim], "unbinding", dimer, MSsoln, NULL , MSnone, 2,(const char**) unbindProducts, outputStates, paramsDe.unbindK);
		}

		smolAddCommandFromString(Reps.sims[jSim], "e ifincmpt all = 0 roiComp stop");
		smolUpdateSim(Reps.sims[jSim]);
		//Reps.sims[jSim]->logfile = nulldev;
	}
}


void dynamicsEngine(simptr currentSim){
	smolRunSimUntil(currentSim, currentSim->time + paramsWe.tau*paramsDe.dt);
}

int findBin(simptr currentSim){
	//Finds number of molecules inside the ROI by inspecting each molecule's location
	int nInBin, nMol;
	double molX, molY;
	nInBin = 0;
	//for(nMolList = 0; nMolList < maxlist;++){ This loop will be included when I add more molecule types to simulations
	for(nMol = 0; nMol < paramsDe.nPart; nMol++){
		molX = currentSim->mols->live[0][nMol]->pos[0]; //first index changes based on how we org mol lists
		molY = currentSim->mols->live[0][nMol]->pos[1];
		if((molX*molX+ molY*molY) < paramsDe.roiR*paramsDe.roiR){
			nInBin++;
		}
	}
	return nInBin;
}

void freeAllSims(){
	int jSim
	for(jSim = 0; jSim <=Reps.iSimMax; jSim++){
	smolFreeSim(Reps.sims[jSim]);
	}
	return 0
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
	double roiParams[] = {0.0, 0.0, 1, 30};
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
	
	Reps.sims[simOut] = smolNewSim(2, lowBounds, highBounds);
	smolSetRandomSeed(Reps.sims[simOut],genrand_int31());
	smolSetGraphicsParams(Reps.sims[simOut], "none", 1, 0);
	smolSetSimTimes(Reps.sims[simOut], 0, 10000, paramsDe.dt); //Reps.sims[simIn]->time is another option
	smolAddMolList(Reps.sims[simOut],"mList");
	smolAddSpecies(Reps.sims[simOut], "A", "mList");
	smolSetSpeciesMobility(Reps.sims[simOut], "A", MSall, paramsDe.difM, NULL, NULL);
	smolSetMaxMolecules(Reps.sims[simOut],paramsDe.nPart);
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
	
	if(paramsDe.reactBit > 0){
			smolAddMolList(Reps.sims[simOut],"dList");
			smolAddSpecies(Reps.sims[simOut], dimer, "dList");
			smolSetSpeciesMobility(Reps.sims[simOut],dimer, MSall, paramsDe.difD, NULL, NULL);
			smolAddReaction(Reps.sims[simOut], "binding", monomer, MSsoln, monomer, MSsoln, 1, dimerptr, &outputStates[0], -1);
			smolSetReactionRate(Reps.sims[simOut], "binding", paramsDe.bindR, 1); //This line allows us to set the rate by the binding radius rather than kOn
			smolAddReaction(Reps.sims[simOut], "unbinding", dimer, MSsoln, NULL , MSnone, 2,(const char**) unbindProducts, outputStates, paramsDe.unbindK);
	}
		
	smolAddCommandFromString(Reps.sims[simOut], "e ifincmpt all = 0 roiComp stop");
	//The following for loop is the major difference between this and the initial configuration simulation building.
	
	
	for(nList = 0; nList < Reps.sims[simIn]->mols->nlist; nList++){
		for(nMol = 0; nMol < Reps.sims[simIn]->mols->nl[nList];nMol++){
			smolAddSolutionMolecules(Reps.sims[simOut], stateList[nList], 1, Reps.sims[simIn]->mols->live[nList][nMol]->pos,Reps.sims[simIn]->mols->live[nList][nMol]->pos);
	}
	}
	smolUpdateSim(Reps.sims[simOut]);
	//Print statement to check functionality. Note molecule indices get changed but everything else preserved
	for(nList = 0; nList < Reps.sims[simIn]->mols->nlist; nList++){
		for(nMol = 0; nMol < Reps.sims[simIn]->mols->nl[nList];nMol++){
		printf("Molecule %i Coords: %f %f \n",nMol, Reps.sims[simOut]->mols->live[nList][nMol]->pos[0],Reps.sims[simOut]->mols->live[nList][nMol]->pos[1]);
		printf("Molecule %i Coords: %f %f \n",nMol, Reps.sims[simIn]->mols->live[nList][nMol]->pos[0],Reps.sims[simIn]->mols->live[nList][nMol]->pos[1]);
		}
	}
}

/*
void copySim1(int simIn, int simOut){

	
	//Takes simIn, reads molecule locations from the sim, then makes a new sim with molecules at the same locations. Most parameters / lines of code are identical to the sim made through initialDist
		
	int nMol;
	//FILE *nulldev = fopen(NULLDEVICE, "w");
	smolFreeSim(Reps.sims[simOut]);
	simSetup(simOut);
	//The following for loop is the major difference between this and the initial configuration simulation building.
	for(nMol = 0; nMol < paramsDe.nPart; nMol++){
		smolAddSolutionMolecules(Reps.sims[simOut], "A", 1, Reps.sims[simIn]->mols->live[0][nMol]->pos,Reps.sims[simIn]->mols->live[0][nMol]->pos);
	}
	smolUpdateSim(Reps.sims[simOut]);
}*/