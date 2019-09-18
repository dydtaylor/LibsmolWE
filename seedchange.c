#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include "twister.c"

/*Seed each simulation based off of rand
Use last ISEED as a reference for possibly no reason
Goal of the program is to be used with the script "prepareRepeatedSim.sh"

Refers to ISEED to initialize C's random number generator.
Iterates the random number generator N times (user arg)
then uses that to reseed ISEED

The goal would be to use each possible "rand" instance at line 49/51 to make a new ISEED for each simulation, but because of how the files are copied with a bash script and each rand call in the loop is "lost", after I call this once I reset the ISEED file once I've copied the new file over, then re-seed with an integer one greater than the previous one.

The reason I'm doing it like this rather than just leaving ISEED the same each time is to hopefully make the runs slightly more reproducible / easier to track the ISEEDs
*/
int main(int argc, char *argv[]){
	int i, j, reset, repeats;
	long iseed, lastiseed;
	FILE *ip, *ip2;
	
	reset=atoi(argv[2]);
	if(reset == 1){
		ip2 = fopen("last ISEED","r");
		fscanf(ip2, "%lx\n", &lastiseed);
		fclose(ip2);
		ip = fopen("ISEED","w");
		fprintf(ip,"%lx\n",lastiseed);
		fclose(ip);
	}
	else{
	ip = fopen("ISEED","r");
	fscanf(ip,"%lx", &iseed);
	fclose(ip);
	
	ip2 = fopen("last ISEED", "w");
	fprintf(ip2, "%lx\n", iseed);
	fclose(ip2);
	
	if(iseed>0) iseed*=-1;
	
	srand(iseed); //Loads ISEED to seed C's rand function
	repeats = atoi(argv[1]);
	for(i = 0; i<repeats-1; i++){
		rand();
	}
	j = rand();
	init_genrand(j);
	
	iseed = genrand_int32();
	ip = fopen("ISEED", "w");
	fprintf(ip, "%lx\n",iseed);
	fclose(ip);
	}
	return 0;
}