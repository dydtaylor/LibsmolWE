9/21/2021
This will provide a basic guide to getting LibsmolWE running and producing data that was used in "Stochastic rare event simulation demonstrates receptor triggeringby kinetic segregation is influenced by oligomerization andclose-contact mechanics" (Taylor, Allard, Read 2021)

# Smoldyn
First, the user will need to install Smoldyn from http://www.smoldyn.org/download.html .
	
* You will need the complete distribution with source code, not the pre-compiled software.
	
* The version of Smoldyn used in the paper is Smoldyn 2.62. Major changes to Smoldyn that might impact these results would be how Smoldyn handles 2D molecular reactions, specifically A+A <-> B reactions. Steve Andrews mentioned in Chapter 2 of his SmolEmulate document, that reaction rate functions have "several problems, some dating back to [the] original 2004 algorithms and some in new algorithms".
	
* The authors of LibsmolWE attempted to select reaction parameters that lie outside of the problematic regions mentioned in "SmolEmulateDoc.pdf", so hopefully minor changes have no impact on the results for 2d binding.
	
* Installation should install the libsmoldyn, smoldyn, smoldynconfigre, libsmoldyn_static.a and libsmoldyn_shared.dylib libraries to the users local folder so that Smoldyn can be called as a library.
	
# Compiling

Update the makefile with your own Smoldyn library locations. The makefile uses dynamic linking rather than static. This should create the executable "weSmoldyn". 

Before running weSmoldyn, set the LD_LIBRARY_PATH variable to the directory where libsmoldyn_static.a is located. e.g.

	LD_LIBRARY_PATH=/YOUR/LIBRARY/LOCATION/HERE/ 
	export LD_LIBRARY_PATH

# Simulation Parameters

 Update the parameters and bin definitions files. Specifically, these files are "binParams.txt", "corralsParams.txt", "dynamicsParams.txt", "WEParams.txt, "binDefinitions.txt".
	
Brief explanation for each line in the parameters
	
1. binDefinitions.txt: Contains the boundaries for each bin, given a WE order parameter of "number of molecules inside the ROI". The first number gives the left-most boundary for the first listed bin, the second number gives the upper boundary for the first listed bin / the lower boundary for the second listed bin, the third gives the upper boundary for the second listed bin / lower boundary for third listed bin, etc. WE is fairly robust to whichever bin definitions you intend to use, though definitions that allow for more replicas nearby the transient bins will tend to allow quicker convergence of flux values. A sample definition for 13 bins is included, and those definitions were used for all simulations where nPart = 256. Note the wide gap for the final bin, which is furthest from the transient bins.
		
1. binParams.txt
			
	1. custom: Whether or not the user intends to use the custom bins in "binDefinitions.txt". If this is 0 then these files are ignored. Publication value was 1.
			
	1. binDims: Number of dimensions for the order parameters. Currently this should only be set to binDims = 1.
			
	1. nBins: Clarifies the number of bins to use in the custom bins. This should be equal to n-1, where n is the number of lines of binDefinitions.txt
		
1. corralsParams.txt: Corrals are currently in an alpha version and as such that functionality is not currently used to produce data in any publications. Results from this functionality have NOT been confirmed through any brute force analysis / analytic analysis / other methods.
			
	1. corralsBit: 1 enables corrals functionality, 0 disables it.
			
	1. corralsWidth: Specifies the width of the corrals created. Corrals are squares tiled, with the center corral being a square with width corralsWidth centered at the origin.
			
	1. corralsRate: Specifies the probabilistic rate that molecules are allowed to pass through the corrals with.
		
1. dynamicsParams.txt
			
	1. dt: Timestep for Smoldyn dynamics simulations. dt = 0.000001 (10^-6) was used. Units of characteristic timescale of system (t*)
			
	1. worldL: Length of square domain, units of the characteristic length scale (L*) of system. Publication values ranged from 3 to 5.333.
			
	1. roiR: Radius of ROI, units of the characteristic length scale of the system (L*). Publication value was always 1.
			
	1. difM: Monomer diffusion coefficient. Units of L*^2/t*. Publication value was always 1.
			
	1. difD: Dimer diffusion coefficient. Units of L*^2/t*. Publication value was always 0.5.
			
	1. bindR: Binding radius of monomers. Units of L*. Publication value was 0.001 (10^-3)
			
	1. unbindK: Unbinding rate of dimers. Units of 1/t*. Publication values ranged from 10^-3 to 10^6.
			
	1. nPart: Maximum number of monomers in the system. Publication values ranged from 2 to 256.
			
	1. reactBit: Bit that enables (1) or disables (0) chemical reactions between monomers.
			
	1. entryBit: Bit that enables (1) or disables (0) probabilistic reentry from close contacts.
			
	1. entryRate: Probabilistic rate that molecules are allowed to enter into the ROI from outside the ROI. Publication values ranged from 0 to 1.

	1. densityBit: When this value is > 0, the worldL given above will be ignored and recalculated based off of nPart (line 8) and density (line 13)
			
	1. density: When densityBit is > 0, this value is used to calculate a new world L based off of nPart. worldL = sqrt(nPart/density). Units of number of molecules per L*^2.
			
	1. monomerStart: If react bit is 1, replicas will either be initialized with purely monomeric (1) solutions with nPart monomers or purely dimeric (0) solutions with nPart/2 dimers.
		
1. WEParams.txt
			
	1. tau: Gives the integer number of Smoldyn timesteps (dt in dynamicsParams.txt) to execute in between WE flux measurement / splitting + merging. Publication value was 50. Units of dt.
			
	1. repsPerBin: AKA mtarg. The number of replicas to maintain in each bin through the splitting / merging process. Publication values ranged from 100 to 200.
			
	1. initialReps: Number of unique replicas to initialize the WE simulation with. These replicas are drawn from a distribution where each molecule is uniformly distributed throughout the entire domain.
			
	1. tauMax: Maximum number of WE steps to use. Only used if FIXEDTIME is defined as 1 in weSmoldyn.h, otherwise simulations run until either the KS tests are passed (see KSTest in weSmoldyn.c) or 250000 seconds of combined computation time has been used. Publication value was 12000, though FIXEDTIME was 0 in all sims for publication. EVEN WITH FIXEDTIME = 0, THIS PARAMETER IS STILL USED: tauMax defines some "default burn in" times that the system will wait for before beginning measurements. For example, when measuring monomerization fractions, the "mCountsWeighted" constitutes a running average that starts after 1/2 of this time has passed.
			
	1. nBins: Primarily useful if there are no custom bin definitions. When there aren't custom bin definitions, this should be equal to nPart in dynamicsParams.txt.
			
	1. fluxBin: Specifies which bin should be used for measuring flux into. This was always 0 in publication.
			
	1. ksNT: How many WE steps should pass before the first KS test is used. This value was always 300 in publication.
			
	1. replaceFluxSims: If (1), then a replica that enters the flux bin will be replaced with a new replica initialized from the starting distribution (uniform molecule distribution). If (0), then a replica that enters the flux bin will have its weight proportionally redistributed to replicas that did not enter the flux bin. This value was 0 for all publication simulations EXCEPT for the close contacts (entryBit = 1 in dynamicsParams.txt), in which case the value was 1.

# Execution

After updating parameters and exporting LD_LIBRARY_PATH, LibsmolWE can now be run using the following command line arguments:
			
1. argv1: ending simfile (Somewhat of a misnomer, gives weight distribution across bins at given timesteps)
			
1. argv2: flux file. Appended to every WE timestep with the flux measured in that timestep
			
1. argv3: Seed/error file. Outputs the ISEED used + miscellaneous errors that I wanted to output here.
			
1. argv4: Save / replace RNG bit. If this is 1 then the RNG seed is saved / ISEED does not increment. If 0	then it increments according to twister.c
			
1. argv5: Execution time file. Outputs, in order: 
	1. Time spent creating initial distribution 
	1. All time spent in the WE splitting/merging 
	1. All time spent in Smoldyn dynamics 
	1. Total time of WE/Dyanamics 
	1. Time spent creating the savestate

1. argv6: Load sim bit. If 1, then savestate.txt is referenced to load the previous savestate, if 0 then this initializes a new WE sim.
					
					
WE Smoldyn outputs the following files:
	
	"Sim File" aka filename given as argv1. This is a cursory overview of how the weight in the WE sim evolves. Output is in series of lines for WE step (tau), then it lists the bin number and the weight within that bin for each bin used. The first 1000 tau steps are always included. After the first 1000 tau steps, it makes a new entry after 0.1% of the current "maximum tau step" has elapsed. Because the maximum tau step can change when the #FIXEDTIME macro is 0, the increments between steps increases the longer the sim runs.
			Currently, the sim file does not output any information not included by "timeSeries.txt"
	
	"Flux file" A line-by-line time series of the flux measured at that given WE step. Each line contains the summed flux measured from each replica that evacuated during that steps.
	
	"Seed file" Outputs the RNG seed (ISEED) used to generate the data. Clarifies what the flux bin was for the data.
	
	"Time file" Line by line, outputs (in seconds) 1. Time spent creating initial distribution 2. All time spent in the WE splitting/merging 3. All time spent in Smoldyn dynamics 4. Total time of WE/Dyanamics 5. Time spent creating the savestate
	
	"savestate.txt": Gives a complete description of the entire WE system at the ending WE step. First line, iSimMax, tells you how many replicas exist in the system. Then, for each replica, it lists: the weight associated with that replica ("Weight w" where "w" is a floating point value), the number of monomers in that replica ("Monomers n" where "n" is an integer), then ordered pairs of all the monomer locations ("X, Y" where x and y are the floating point coordinates in 2d). Finally, it lists how many dimers are in the system ("Dimers m" where "m" is an integer) followed by ordered pairs of all the dimer locations ("X, Y" where x and y are the floating point coordinates in 2d).
	
	"timeSeries.txt": Each line is in the following format: "B# W mFrac tau", B# is the bin number, W is the total weight inside that bin, mFrac is the weight-adjusted monomerization fraction of replicas inside that bin, and tau is the WE step that the 3 prior numbers were measured at.
	
	"mCountsWeighted.txt": The first line gives the weighted number of monomers measured (SUM(weight_i * nMonomers_i)), the number of dimers measured (SUM(weight_i * nDimers_i)), and the total weight measured (SUM(weight_i)) beginning at a timestep tauMax/2. The sequential lines give a running sum of the previous line with the weighted number of monomers, weighted number of dimers, and total weight measured at that time step. To calculate the monomerization fraction from one of the lines, use nMonomers_weighted / (nPart * total_weight), i.e. (First entry of the line) / (3rd entry of the same line * nPart)
	
	ksOut.txt: Gives the calculations used to measure the KS statistic of a run. This should consist of 5 line blocks listing the maximum flux measured to create the KS statistic, the empirical histogram for the middle third of flux data, the empirical histogram for the final third of flux data, the KS statistic measured from these histograms, and the number of WE steps used for the histogram (nT). a ksStat of 2 means that the data did not have enough non-0 measurements. (line 33, #NNZMIN in weSmoldyn.h)
	
	dualKS.txt: Same as ksOut.txt except it modifies the data used to create the histograms slightly. It removes a number of "0s" from the 0 histogram bin equal to the minimum number of zeroes measured in either third. E.g. if the middle third of data has 1000 measurements and 500 of them are 0s, and the final third of data has 1000 measurements and 300 of them are 0s, then the file gives a KSstat comparing 700 measurements from the middle third and 700 measurements from the final third, where the 300 measurements that are omitted are all 0s. There is an additional line included below nT that gives the number of zeroes counted from the middle and final third as well as confirmation of the number of measurements removed to make the dksStat.
	
	The files included in this distribution will allow you to execute weSmoldyn without changing any of the parameters to obtain data for L = 5.333 and N = 256 monomer-only simulations. The terminal command is as follows:
	
	./weSmoldyn sampleOut.txt sampleFlux.txt sampleSeed.txt 0 sampletime.txt 0
	
	Smoldyn outputs a lot of junk to stdout that is not useful outside of debugging purposes. For data collection, I would recommend suppressing that output by including &>/dev/null at the end of the terminal command. However, verifying that Smoldyn is doing what you want it to is useful, in which case I would recommend using a debugger like lldb and creating a breakpoint at smolDynamics.c:123 (alternatively, smolDynamics.c:564 if you want one copied during the splitting), where the smolDisplaySim command is. This will provide all details of the specific replica initialized, which should be identical to the other initialized replicas. Debugging Smoldyn features in LibsmolWE is cumbersome, so I would recommend using this smolDisplaySim output as a comparison tool in combination with vanilla Smoldyn / Libsmol simulations.
	
	To calculate the evacuation time / MFPT from the flux, the formula is given by the Hill relation, MFPT = 1 / mean(flux). It is common practice to discard the first few (1/4 to 1/2 of the data, though discarding less may be fine) flux measurements to account for shifts from the initial distribution to the non-equilibrium steady-state distribution of weights obtained during the WE simulation.