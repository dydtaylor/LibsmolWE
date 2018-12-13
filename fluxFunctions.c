void newSim(int iSim){

	int nMol;
	
	double lowBounds[] = {-paramsDe.worldLength/2,-paramsDe.worldLength/2};
	double highBounds[] = {paramsDe.worldLength/2, paramsDe.worldLength/2};
	double botLeftCornerRect[] = {-paramsDe.worldLength/2, -paramsDe.worldLength/2, paramsDe.worldLength};
	double topRightCornerRect[] = {paramsDe.worldLength/2, paramsDe.worldLength/2, -paramsDe.worldLength};
	double roiParams[] = {0.0, 0.0, 1, 30};
	double insideRoi[] = {0.0, 0.0};
	
	Reps.sims[iSim] = smolNewSim(2, lowBounds, highBounds);
	smolSetRandomSeed(Reps.sims[iSim],genrand_int31());
	smolSetGraphicsParams(Reps.sims[iSim], "none", 1, 0);
	smolSetSimTimes(Reps.sims[iSim], 0, 10000, paramsDe.dt); //Reps.sims[simIn]->time is another option
	smolAddSpecies(Reps.sims[iSim], "A", NULL);
	smolSetSpeciesMobility(Reps.sims[iSim], "A", MSall, paramsDe.difC,0,0);
	smolSetMaxMolecules(Reps.sims[iSim],paramsDe.nPart);
	smolAddSurface(Reps.sims[iSim],"bounds");
	smolAddPanel(Reps.sims[iSim], "bounds", PSrect, NULL, "-x", topRightCornerRect);
	smolAddPanel(Reps.sims[iSim], "bounds", PSrect, NULL, "-y", botLeftCornerRect);
	smolAddPanel(Reps.sims[iSim], "bounds", PSrect, NULL, "+x", botLeftCornerRect);
	smolAddPanel(Reps.sims[iSim], "bounds", PSrect, NULL, "+y", topRightCornerRect);
	smolSetSurfaceAction(Reps.sims[iSim], "bounds", PFboth, "all", MSall, SAreflect);
		
	//ROI Surface + compartment
	smolAddSurface(Reps.sims[iSim], "roi");
	smolAddPanel(Reps.sims[iSim],"roi", PSsph, NULL, "+0", roiParams);
	smolSetSurfaceAction(Reps.sims[iSim], "roi", PFboth, "all", MSall, SAtrans);
	smolAddCompartment(Reps.sims[iSim],"roiComp");
	smolAddCompartmentSurface(Reps.sims[iSim],"roiComp","roi");
	smolAddCompartmentPoint(Reps.sims[iSim],"roiComp",insideRoi);
	smolAddSolutionMolecules(Reps.sims[iSim],"A",paramsDe.nPart, lowBounds, highBounds);
}

double fluxNewSim(){
	double fluxOut = 0;
	int nFlux = Reps.binContentsMax[paramsWe.fluxBin];
	if (nFlux>0){
		for(iReps = nFlux -1; iReps >= 0; iReps--){
			fluxOut += Reps.weights[Reps.binContents[iReps][paramsWe.fluxBin]];
			smolFreeSim(Reps.sims[Reps.binContents[iReps][paramsWe.fluxBin]]);
			
		}
	}
}