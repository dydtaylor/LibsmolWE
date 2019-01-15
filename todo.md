## libsmolWE code review request list

### Folder structure

* Make a folder called /analysis and put all the m-files into it.

* Make a folder called /scripts and put all the .sh scripts into it.

* Make a folder __outside the git repo__ called /runs. This is where simulations should happen, and where any saved output from analysis should go.

* Add a README.md (markdown) with the basic command-line call. I have started this in my pull-request. Add to my outline.


### Algorithm and code

* Something seems funny in flux(), for the replica replacement. After a flux event, do you copy the replica with the largest iSim? Is that the correct algorithm? I think you just do nothing, and wait for the next splitting cycle to create a new replica. Or else you pick a replica at random, weighted by replica weight, and copy that one. I forget which one it is supposed to be.

* What is the array (0,0,1,30) in smolSimCopy and smolDynamics? Should the number 30 be a macro rather than hard-coded in?

* Use a macro for 10000, the smoldyn maximum time per replica (I think).

* Move simCopy1() into smolDynamics.c. Create a function called simSetup() or similar, and put L1-36 (which is identical to a block in simCopy()) into it. Make the parameters global by putting them at the top of smolDynamics.c, outside the function. Then, call simSetup() inside initialDist() and simCopy1(), instead of duplicating code.

* Given how similar the first-quarter and last-3-quarters are, consider just putting the flux() into an if statement that checks whether you're in the first quarter or not. Any downsides?


* Suppress stdout, if possible.

* Move parameter input to a new function getParams() (which can still be in weSmoldyn.c, just not in main())

* Make parameter input files more human-readable: Put each parameter on a separate line, and make C code ignore all text after the first word (the numerical value). Then in the txt file, you can put a description of the parameter in text, which will be ignored by C.

* there are many empty if statements that appear to be debugging assertions (sanity checks). Put these in a statement with DEBUGGING macro
```
#DEFINE DEBUGGING 1
if(DEBUGGING && your condition)
```

### Goals mentioned in quarter report

* Add KS test into code

* Generate MFPTs and plot against brute force
