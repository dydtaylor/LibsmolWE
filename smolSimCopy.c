void copySim1(int simIn, int simOut){

	int nMol;
	//FILE *nulldev = fopen(NULLDEVICE, "w");
	
	double lowBounds[] = {-paramsDe.worldLength/2,-paramsDe.worldLength/2};
	double highBounds[] = {paramsDe.worldLength/2, paramsDe.worldLength/2};
	double botLeftCornerRect[] = {-paramsDe.worldLength/2, -paramsDe.worldLength/2, paramsDe.worldLength};
	double topRightCornerRect[] = {paramsDe.worldLength/2, paramsDe.worldLength/2, -paramsDe.worldLength};
	double roiParams[] = {0.0, 0.0, 1, 30};
	double insideRoi[] = {0.0, 0.0};
	
	Reps.sims[simOut] = smolNewSim(2, lowBounds, highBounds);
	smolSetRandomSeed(Reps.sims[simOut],genrand_int31());
	smolSetGraphicsParams(Reps.sims[simOut], "none", 1, 0);
	smolSetSimTimes(Reps.sims[simOut], 0, 10000, paramsDe.dt); //Reps.sims[simIn]->time is another option
	smolAddSpecies(Reps.sims[simOut], "A", NULL);
	smolSetSpeciesMobility(Reps.sims[simOut], "A", MSall, paramsDe.difC,0,0);
	smolSetMaxMolecules(Reps.sims[simOut],paramsDe.nPart);
	smolAddSurface(Reps.sims[simOut],"bounds");
	smolAddPanel(Reps.sims[simOut], "bounds", PSrect, NULL, "-x", topRightCornerRect);
	smolAddPanel(Reps.sims[simOut], "bounds", PSrect, NULL, "-y", botLeftCornerRect);
	smolAddPanel(Reps.sims[simOut], "bounds", PSrect, NULL, "+x", botLeftCornerRect);
	smolAddPanel(Reps.sims[simOut], "bounds", PSrect, NULL, "+y", topRightCornerRect);
	smolSetSurfaceAction(Reps.sims[simOut], "bounds", PFboth, "all", MSall, SAreflect);
		
	//ROI Surface + compartment
	smolAddSurface(Reps.sims[simOut], "roi");
	smolAddPanel(Reps.sims[simOut],"roi", PSsph, NULL, "+0", roiParams);
	smolSetSurfaceAction(Reps.sims[simOut], "roi", PFboth, "all", MSall, SAtrans);
	smolAddCompartment(Reps.sims[simOut],"roiComp");
	smolAddCompartmentSurface(Reps.sims[simOut],"roiComp","roi");
	smolAddCompartmentPoint(Reps.sims[simOut],"roiComp",insideRoi);
		
	smolAddCommandFromString(Reps.sims[simOut], "e ifincmpt A = 0 roiComp stop");
	for(nMol = 0; nMol < paramsDe.nPart; nMol++){
		smolAddSolutionMolecules(Reps.sims[simOut], "A", 1, Reps.sims[simIn]->mols->live[0][nMol]->pos,Reps.sims[simIn]->mols->live[0][nMol]->pos);
	}
	smolUpdateSim(Reps.sims[simOut]);
	//Reps.sims[simOut]->logfile = nulldev;
	for(nMol = 0; nMol < paramsDe.nPart; nMol++){
		//printf("Molecule %i Coords: %f %f \n",nMol, Reps.sims[simOut]->mols->live[0][nMol]->pos[0],Reps.sims[simOut]->mols->live[0][nMol]->pos[1]);
		//printf("Molecule %i Coords: %f %f \n",nMol, Reps.sims[simIn]->mols->live[0][nMol]->pos[0],Reps.sims[simIn]->mols->live[0][nMol]->pos[1]);
	}
}
/*
void copySim2(simptr simIn, simptr simOut){
	FILE *simFile;
	simFile = fopen("simStore.txt","w");
	writesim(simIn, simFile);
	fclose(simFile);
	simOut = smolPrepareSimFromFile(NULL, "simStore.txt", NULL);
}*/