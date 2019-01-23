void simSetup(int simOut){

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

void copySim1(int simIn, int simOut){

	/*
	Takes simIn, reads molecule locations from the sim, then makes a new sim with molecules at the same locations. Most parameters / lines of code are identical to the sim made through initialDist
	*/
	
	int nMol;
	//FILE *nulldev = fopen(NULLDEVICE, "w");
	smolFreeSim(Reps.sims[simOut]);
	simSetup(simOut);
	//The following for loop is the major difference between this and the initial configuration simulation building.
	for(nMol = 0; nMol < paramsDe.nPart; nMol++){
		smolAddSolutionMolecules(Reps.sims[simOut], "A", 1, Reps.sims[simIn]->mols->live[0][nMol]->pos,Reps.sims[simIn]->mols->live[0][nMol]->pos);
	}
	smolUpdateSim(Reps.sims[simOut]);
}