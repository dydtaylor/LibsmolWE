/**** Compiling and linking ****

Compile source code to object code with:
gcc -Wall -O0 -g -c testcode.c

Static link, if Libsmoldyn was compiled without OpenGL:
gcc testcode.o /usr/local/lib/libsmoldyn_static.a -o testcode

Dynamic link, if Libsmoldyn was compiled without OpenGL:
gcc testcode.o -o testcode -lsmoldyn_shared

Static link, if Libsmoldyn was compiled with OpenGL:
gcc testcode.o /usr/local/lib/libsmoldyn_static.a -I/System/Library/Frameworks/OpenGL.framework/Headers -I/System/Library/Frameworks/GLUT.framework/Headers -framework GLUT -framework OpenGL -framework Cocoa -L/System/Library/Frameworks/OpenGL.framework/Libraries -o testcode -ltiff

Dynamic link, if Libsmoldyn was compiled with OpenGL:
gcc testcode.o -I/System/Library/Frameworks/OpenGL.framework/Headers -I/System/Library/Frameworks/GLUT.framework/Headers -framework GLUT -framework OpenGL -framework Cocoa -L/System/Library/Frameworks/OpenGL.framework/Libraries -o testcode -lsmoldyn_shared -ltiff
*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "libsmoldyn.h"			// This line is required for accessing Libsmoldyn
#include <vector>
#include <fstream>
#include <iterator>
#include <string>
#include <iostream>
#include <time.h>
#include "mex.h"

simptr getNewSim(simptr sim,double* posParticles, double* species, int seed,double tau,double kon , double koff, double difcm, double difcd,double dt,double stopflag, double radiusROI, double boxLength)
{
// returns a pointer to a new simulation initialized from the input parameters  

//  if (sim==NULL)
 // {
  //  smolFreeSim(sim);
  //} 
  

	double v1[2],v2[2],color[4],t;
	char **products;
	int ctr,ctf;

	//sim=NULL;
	//printf("Test program.\n");

	//products=(char**) calloc(2,sizeof(char*));
	//products[0]=(char*) calloc(256,sizeof(char));
	//products[1]=(char*) calloc(256,sizeof(char));

	//smolSetDebugMode(1);
	v1[0]=-boxLength/2;
	v1[1]=-boxLength/2;
	v2[0]=boxLength/2;
	v2[1]=boxLength/2;
	sim=smolNewSim(2,v1,v2);
	smolSetSimTimes(sim,0,tau,dt);
	smolSetBoundaryType(sim,-1,-1,'r');
	smolAddMolList(sim,"A");
	smolAddMolList(sim,"AA");
	smolAddSpecies(sim,"A","A");
    smolAddSpecies(sim,"AA","AA");

	//smolAddSpecies(sim,"fox","flist");
	smolSetSpeciesMobility(sim,"A",MSall,difcm,NULL,NULL);
    smolSetSpeciesMobility(sim,"AA",MSall,difcd,NULL,NULL);
    const char *dimerization = "dimerization";
    const char *undimerization = "undimerization";
    const char *A = "A";
    const char *AA = "AA";
    enum MolecState stateA ;
    stateA= MSsoln; 
    enum MolecState stateAA ;
    stateAA = MSsoln; 

    int nproduct = 1;
    char **prodspeciesA;
    char **prodspeciesAA;
    prodspeciesA=(char**) calloc(2,sizeof(char*));
	prodspeciesA[0]=(char*) calloc(256,sizeof(char));
	prodspeciesA[1]=(char*) calloc(256,sizeof(char));
    prodspeciesAA=(char**) calloc(1,sizeof(char*));
	prodspeciesAA[0]=(char*) calloc(256,sizeof(char));



    strcpy(prodspeciesA[0],"A");
	strcpy(prodspeciesA[1],"A");
    strcpy(prodspeciesAA[0],"AA");




	enum MolecState prodstatesAA[1];
    prodstatesAA[0] = MSsoln;
    enum MolecState prodstatesA[2];
    prodstatesA[0] = MSsoln;
    prodstatesA[1] = MSsoln;





    smolAddReaction(sim,"dimerization","A",stateA,"A",stateA,1,(const char**) prodspeciesAA, prodstatesAA,kon);
    smolAddReaction(sim,"undimerization","AA",stateAA,NULL,MSnone,2,(const char**) prodspeciesA, prodstatesA,koff);


	//prodstates[0]=MSsoln;
	//prodstates[1]=MSsoln;
	//smolAddReaction(sim,"r1","rabbit",MSsoln,NULL,MSnone,2,(const char**) products,prodstates,10);
	//strcpy(products[0],"fox");
	//strcpy(products[1],"fox");
	//smolAddReaction(sim,"r2","rabbit",MSsoln,"fox",MSsoln,2,(const char**) products,prodstates,8000);
	//smolAddReaction(sim,"r3","fox",MSsoln,NULL,MSnone,0,NULL,NULL,10);

  for (int p = 0; p<( (int) species[0]*2); p+=2)
  {
    double v3[] = {posParticles[p],posParticles[p+1]};
	  smolAddSolutionMolecules(sim,"A",1,v3,v3);
  }
  for (int p = (int)  species[0]*2; p<(( (int) species[0]*2)+( (int) species[1]*2)); p+=2)
  {
    double v3[] = {posParticles[p],posParticles[p+1]};
	  smolAddSolutionMolecules(sim,"AA",1,v3,v3);
  }
	//smolAddSolutionMolecules(sim,"fox",1000,v1,v2);


    // compartment box

    smolAddSurface(sim, "boxwalls");    
    double p1[] = {-boxLength/2, -boxLength/2, boxLength};
    double p2[] = {boxLength/2, boxLength/2, -boxLength};
    double p3[] = {-boxLength/2, boxLength/2, boxLength};
    double p4[] = {boxLength/2, -boxLength/2, -boxLength};



    smolAddPanel(sim, "boxwalls", PSrect, NULL , "+x" , p1);
    smolAddPanel(sim, "boxwalls", PSrect, NULL , "-x" , p2);
    smolAddPanel(sim, "boxwalls", PSrect, NULL , "+y" , p3);
    smolAddPanel(sim, "boxwalls", PSrect, NULL , "-y" , p4);

    smolSetSurfaceAction(sim, "boxwalls", PFboth, "all", MSall, SAreflect);


    smolAddCompartment(sim,"boxcompartment");
    smolAddCompartmentSurface(sim,"boxcompartment", "boxwalls");
    double p5[] = {0.0,0.0};
    smolAddCompartmentPoint(sim,"boxcompartment",p5);


    // region of interest

    smolAddSurface(sim, "roi");    
    double p6[] = {0.0, 0.0, radiusROI, radiusROI/2};



    smolAddPanel(sim, "roi", PSsph, NULL , "+0" , p6);

    smolSetSurfaceAction(sim, "roi", PFboth, "all", MSall, SAtrans);


    smolAddCompartment(sim,"roicompartment");
    smolAddCompartmentSurface(sim,"roicompartment", "roi");
    double p7[] = {0.0,0.0};
    smolAddCompartmentPoint(sim,"roicompartment",p7);
    int m,ll;
    // stop conditon
    //smolSetOutputPath(sim,"./");
    //smolAddOutputFile(sim,"F", -1, 1);
    //smolAddCommandFromString(sim,"e ifless A 2 stop");
    if (stopflag >0)
    {
        smolAddCommandFromString(sim,"e ifincmpt all < 1 roicompartment stop");
    }
    //smolAddCommandFromString(sim,"a molcountincmpt roicompartment FPT.txt");
    //smolAddCommandFromString(sim,"N 10 molcount F");

    //smolAddOutputFile(sim,"F", -1, 1);


    // simulation start
    
    smolSetRandomSeed(sim,seed);


    smolUpdateSim(sim);
    //smolDisplaySim(sim);

    return sim;


}




void DynamicsEngine(double *inMatrix, double *outMatrix , double *outMatrix2 , double *outMatrix3,double  *numParticles,double tau,double kon, double koff, double difcm, double difcd, int seed,double dt,double stopflag, double radiusROI, double boxLength) {
// takes input configuration and simulates for one timestep and returns output 
// outMatrix - trajectory after tau
// outMatrix2 - number of A and AA after tau


  //srand (time (NULL));
  //char **spname,string[STRCHAR];
  int d,dim,ll,m,t;
  simptr sim =NULL ; 

  std::vector<double> times;
     // simptr sim ; 
      //int seed = rand();
      sim = getNewSim(sim,inMatrix,numParticles,seed,tau,kon,koff,difcm,difcd,dt,stopflag,radiusROI,boxLength);

      smolRunSim(sim);
      molssptr mols;
      moleculeptr mptr;
      mols=sim->mols;
      //printf("n=%d",mols->nlist);
      //printf("nl=%d",mols->nl[0]);
      //printf("time=%g\n",sim->time);
      times.push_back(sim->time);
      int outMatrixDim = 0;
      int prev = 0;
      int ind;
      for(ll=0;ll<mols->nlist;ll++){
         //char mollist[100];
         //printf(smolGetMolListName(sim,ll,mollist));
         //printf("\n");
         //printf("size:%d",(int) mols->nl[ll]);
        for(m=0;m<mols->nl[ll];m++) {
          //printf("ll=%d",ll);
          //printf("m=%d\n",m);
          mptr=mols->live[ll][m];
          dim=sim->dim;
          //printf("dim=%d",dim);
          for(d=0;d<dim;d++) {
          //printf("pos= %g\n",mptr->pos[d]);
          double pos = mptr->pos[d];
          outMatrix[d+(m*2)+prev] = pos;
          ind = d+(m*2); 
          //char mollist[100];
          //printf(smolGetMolListName(sim,ll,mollist));
          //printf("\n");
          //printf("ind:%g",(int) ind);
          //printf("r:%g \n",outMatrix[d+(m*2)]);
        }
        

      }
      prev = ind+1;


   } 
   outMatrix2[0] = (double)  mols->nl[0];
   outMatrix2[1] = (double)  mols->nl[1];
   outMatrix3[0] = (double)  sim->time;


    smolFreeSim(sim);






	}




void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

    // a function that converts the matlab function into a C++ one by converting data types and preparing output data types


    double *inMatrix;               /* 1xN input matrix */
    double nrows;                   /* size of matrix */

    double *outMatrix;              /* output matrix */
    double *outMatrix2;              /* output matrix */
    double *outMatrix3;
    double *numParticles;
    double *numberParticles;
    double *difcd;
    double *difcm;
    double *radiusROI; 			/*Region of Interest Radius*/
    double *boxLength;				/*Box Length*/
    int seed;
    double *tau;
    double *kon;
    double *koff;
    double *dt;
    double stopflag;
    difcm = mxGetPr(mxGetField(prhs[3],0,"difcm"));
    difcd = mxGetPr(mxGetField(prhs[3],0,"difcd"));
    tau  = mxGetPr(mxGetField(prhs[2],0,"tau"));
    kon  = mxGetPr(mxGetField(prhs[3],0,"kon"));
    koff  = mxGetPr(mxGetField(prhs[3],0,"koff"));
    dt = mxGetPr(mxGetField(prhs[4],0,"dt"));
    numberParticles = mxGetPr(mxGetField(prhs[4],0,"numberParticles"));
    radiusROI = mxGetPr(mxGetField(prhs[3],0,"radiusROI"));
    boxLength = mxGetPr(mxGetField(prhs[3],0,"boxLength"));
    seed = (int) mxGetScalar(prhs[5]);
    stopflag = mxGetScalar(prhs[5]);

    inMatrix = mxGetPr(prhs[0]);
    numParticles = mxGetPr(prhs[1]);


    printf("NP:%f",  numberParticles[0]);
    printf("tau:%f", tau[0]);
    printf("koff:%f", koff[0]);
    printf("kon:%f", kon[0]);
    printf("dt:%f", dt[0]);
    printf("seed:%d", seed);






    double ncols = numberParticles[0]*2;
  //  int sizem = nrows*ncols;
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)ncols,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,2,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);



    //plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);

    /* get a pointer to the real data in the output matrix */
     outMatrix = mxGetPr(plhs[0]);
     outMatrix2 = mxGetPr(plhs[1]); // A and AA
     outMatrix3 = mxGetPr(plhs[2]);
  // outMatrix[0] =  sizem;  

    

    DynamicsEngine(inMatrix,outMatrix,outMatrix2,outMatrix3,numParticles,tau[0],kon[0],koff[0],difcm[0],difcd[0],seed,dt[0],stopflag,radiusROI[0],boxLength[0]);


    //evacuations(inMatrix,outMatrix,numParticles,(int) sizem, (int) ncols);


}






