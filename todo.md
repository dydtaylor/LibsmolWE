## libmolWE code review request list

* Something seems funny in flux(), for the replica replacement. Do you just copy the replica with the largest iSim? I don't think that's the correct algorithm. I think you just do nothing, and wait for the next splitting cycle to create a new replica.

* What is the array (0,0,1,30) in smolSimCopy and smolDynamics? Should this be a macro rather than hard-coded in?
* Use a macro for 10000, the smoldyn maximum time per replica (I think).


* Add a README.md (markdown) with the basic command-line call

* Suppress stdout, if possible.

* Move parameter input to a new function getParams(), which can still be in weSmoldyn.c, just not in main()
* Make parameter input files more human-readable: Put each parameter on a separate line, and make C code ignore all text after the first word (the numerical value). Then in the txt file, you can put a description of the parameter in text, which will be ignored by C.

* there are many empty if statements that appear to be debugging assertions (sanity checks). Put these in a statement with DEBUGGING macro
#DEFINE DEBUGGING 1
if(DEBUGGING && your condition)




### Goals mentioned in quarter report

* Add KS test
* Generate MFPTs and plot against brute force
