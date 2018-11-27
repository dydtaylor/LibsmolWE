#include "mfpt.h"


void initialDist(int nInit){
	int jSim;
	double lowBounds[] = {-paramsDe.worldLength/2,-paramsDe.worldLength/2};
	double highBounds[] = {paramsDe.worldLength/2, paramsDe.worldLength/2};
	double botLeftCornerRect[] = {-paramsDe.worldLength/2, -paramsDe.worldLength/2, paramsDe.worldLength};
	double topRightCornerRect[] = {paramsDe.worldLength/2, paramsDe.worldLength/2, -paramsDe.worldLength};
	double roiParams[] = {0.0, 0.0, 1, 30};
	double insideRoi[] = {0.0, 0.0};
	for(jSim = 0; jSim < nInit; jSim++){
		//Molecules + BCs
		Reps.sims[jSim] = smolNewSim(2, lowBounds,highBounds);
		smolSetRandomSeed(Reps.sims[jSim],genrand_int31());
		smolSetGraphicsParams(Reps.sims[jSim], "opengl", 1, 0);
		smolSetSimTimes(Reps.sims[jSim],0,10000,paramsDe.dt);
		smolAddSpecies(Reps.sims[jSim],"A",NULL);
		smolSetSpeciesMobility(Reps.sims[jSim],"A",MSall, paramsDe.difC, 0, 0);
		smolAddSolutionMolecules(Reps.sims[jSim], "A", paramsDe.nPart, NULL, NULL);
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
		
		smolAddCommandFromString(Reps.sims[jSim], "e ifincmpt A = 0 roiComp stop");
		smolUpdateSim(Reps.sims[jSim]);
	}
}

void dynamicsEngine(simptr currentSim){
	double currentT, breakT;
	
	currentT = currentSim->time;
	breakT = currentT + paramsWe.tau*paramsDe.dt;
	smolRunSimUntil(currentSim, breakT);
}