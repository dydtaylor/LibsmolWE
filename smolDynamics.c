void initialDist(int nInit){
/*
!---------------------------------------------------------------------------------------------------------------------
!	Description:
!		This function creates the initial distribution of smoldyn replicas. Each replica is created in a monomer-only
!		configuration where the monomers are distributed randomly throughout the system.
!	Method:
!		1.Create several local variables for dynamics parameters, reaction product states, and molecule string labels
!			reformatted for smoldyn functions
!		2.For the empty simptrs in Reps that will be initialized, create new Smoldyn sims. Each Sim is created by:
!			a: Create new sim, set the random seed, disable graphics, set sim times
!			b: Create list to store molecules in, define the molecules, which list they belong in, set their Diffusion
!				and specify the maximum number of molecules. If reactions are enabled, we do this for dimers as well
!				and then add the reaction in and set the reaction rates. Afterwards we place nPart monomers randomly
!			c: Define a surface for boundaries of the system. Specify reflective boundary conditions.
!			d: Define a surface and compartment for the region of interest. Neither directly interacts with the mols
!			e: Specify that simulations should terminate once all molecules evacuate from the RoI
!		
!---------------------------------------------------------------------------------------------------------------------
*/
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
		smolSetSimTimes(Reps.sims[jSim],0,SMOLTIMEMAX,paramsDe.dt);
		
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
		
		
		//ROI Surface + compartment
		smolAddSurface(Reps.sims[jSim], "roi");
		smolAddPanel(Reps.sims[jSim],"roi", PSsph, NULL, "+0", roiParams);
		if(!paramsDe.reentryRateBit){
		smolSetSurfaceAction(Reps.sims[jSim], "roi", PFboth, "all", MSall, SAtrans);
		}
		if(paramsDe.reentryRateBit){
			smolSetSurfaceAction(Reps.sims[jSim], "roi", PFback, "all", MSall, SAtrans);
			smolSetSurfaceAction(Reps.sims[jSim], "roi", PFfront, "all", MSall, SAreflect);
			smolSetSurfaceRate(Reps.sims[jSim], "roi", "all", MSsoln, MSsoln, MSbsoln, paramsDe.reentryRate, NULL, 1);
		}
		smolAddCompartment(Reps.sims[jSim],"roiComp");
		smolAddCompartmentSurface(Reps.sims[jSim],"roiComp","roi");
		smolAddCompartmentPoint(Reps.sims[jSim],"roiComp",insideRoi);
		smolUpdateSim(Reps.sims[jSim]);
		smolDisplaySim(Reps.sims[jSim]);
		if(STOPCOMMAND){
		smolAddCommandFromString(Reps.sims[jSim], "e ifincmpt all = 0 roiComp stop");
		}
		
		//Reps.sims[jSim]->logfile = nulldev;
	}
}


void dynamicsEngine(simptr currentSim){
/*
!---------------------------------------------------------------------------------------------------------------------
!		Description:
!			This function executes all of the dynamics for a single simptr through an entire tau step.
!		Method:
!			Receives a simptr as input. Runs a single smoldyn timestep on that sim for each tau.
!---------------------------------------------------------------------------------------------------------------------
*/
	//smolUpdateSim(currentSim);
	for(int i = 0; i < paramsWe.tau; i++){
	smolRunTimeStep(currentSim);
	}
	
}

int findBin(simptr currentSim){
/*
!---------------------------------------------------------------------------------------------------------------------
!		Description:
!			This examines a simptr then determines the bin it should be in based on the order parameter used
!		Method:
!			If ROBINS macro is enabled, it loops through all molecule lists and increments the bin every time it finds
!			a molecule inside the RoI radius (x^2 + y^2 < R^2). If not using ROBINS, the bin out is the number of live
!			monomers.
!---------------------------------------------------------------------------------------------------------------------
*/
	//Finds number of molecules inside the ROI by inspecting each molecule's location
	int binOut, nMol, nMolList;
	double molX, molY;
	binOut = 0;
	if(ROBINS){
	for(nMolList = 0; nMolList < currentSim->mols->nlist; nMolList++){ 
	for(nMol = 0; nMol < currentSim->mols->nl[nMolList]; nMol++){
		molX = currentSim->mols->live[nMolList][nMol]->pos[0]; //first index changes based on how we org mol lists
		molY = currentSim->mols->live[nMolList][nMol]->pos[1];
		if((molX*molX+ molY*molY) < paramsDe.roiR*paramsDe.roiR){

			binOut++;
		}
	}
	}
	};
	if(!ROBINS){
		binOut = currentSim->mols->nl[0]; //Sets bin to # of monomers
	}
	return binOut;
}

void freeAllSims(){
/*
!---------------------------------------------------------------------------------------------------------------------
!		Description:
!			This frees all smoldyn simulations.
!		Method:
!			Executes a loop of "smolFreeSim" for each sim currently being used.
!---------------------------------------------------------------------------------------------------------------------
*/
	int jSim;
	for(jSim = 0; jSim <=Reps.iSimMax; jSim++){
	smolFreeSim(Reps.sims[jSim]);
	}
}

void copySim1(int simIn, int simOut){

/*
!---------------------------------------------------------------------------------------------------------------------
!		Description:
!			When a sim is to be copied (split) this creates a new sim
!		Method:
!			This functionally is very similar to the function initialDist. There are two key differences:
!			1. Rather than creating multiple sims, this only creates one sim.
!			2. The molecule positions are not placed randomly, but at the same location as a previous replica's. The
!			positions are copied from Reps.sims[simIn] and the output replica is Reps.sims[simOut].
!---------------------------------------------------------------------------------------------------------------------
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
		if(!paramsDe.reentryRateBit){
		smolSetSurfaceAction(Reps.sims[simOut], "roi", PFboth, "all", MSall, SAtrans);
		}
		if(paramsDe.reentryRateBit){
			smolSetSurfaceAction(Reps.sims[simOut], "roi", PFback, "all", MSall, SAtrans);
			smolSetSurfaceAction(Reps.sims[simOut], "roi", PFfront, "all", MSall, SAreflect);
			smolSetSurfaceRate(Reps.sims[simOut], "roi", "all", MSsoln, MSsoln, MSbsoln, paramsDe.reentryRate, NULL, 1);
		}
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

