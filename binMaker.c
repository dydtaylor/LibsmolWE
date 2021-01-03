#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include "twister.c"

#define RAND genrand_real1()
#define RANDINT genrand_int32()
#define PI 3.14159265

int main(int argc, char *argv[]){
	
	/*Command line inputs:
	argv[1]: Number of molecules
	argv[2]: roi Radius
	argv[3]: domain L
	*/
	//Input variables
	int n, nBins, j;
	double L, roiR;
	n = atoi(argv[1]);
	sscanf(argv[2],"%lf",&roiR);
	sscanf(argv[3],"%lf",&L);
	//Calculation variables
	double p, meanIn,sd;
	p = roiR*roiR*PI/(L*L);
	meanIn = (double) p*n;
	sd = sqrt((double)n*(1-p)*p);
	//File variables
	FILE *binFile = fopen("binDefinitions.txt","w");
	FILE *binParams = fopen("binParams.txt","w");
	
	//Set up bin bounds
	double bin3sd, bin2sd, bin1sd;
	bin3sd = meanIn-3*sd;
	bin2sd = meanIn - 2*sd;
	bin1sd = meanIn - sd;
	
	if(bin1sd < 1){
		bin1sd = 1;
	}
	int nEveryMax; //Where the fine-grained bins stop
	nEveryMax = floor(bin3sd);
	bin1sd = ceil(bin1sd);
	bin2sd = ceil(bin2sd);
	bin3sd = ceil(bin3sd);
	meanIn = ceil(meanIn);
	
	nBins = -1;
	for(j = 0; j < nEveryMax + 1; j++){
		fprintf(binFile,"%i\n",j);
		nBins++;
	}
	fprintf(binFile,"%i\n",(int)bin3sd);
	fprintf(binFile,"%i\n",(int)bin2sd);
	fprintf(binFile,"%i\n",(int)bin1sd);
	fprintf(binFile,"%i\n",(int)meanIn);
	fprintf(binFile,"%i\n",n);
	nBins = nBins+5;
	
	fprintf(binParams,"custom\t1\nbinDims\t1\nnBins\t%i",nBins);
	fclose(binFile);
	fclose(binParams);
	return 0;
}