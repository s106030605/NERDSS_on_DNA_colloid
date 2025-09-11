# nerdss
Software for structure-resolved RD

1. To compile the code and install the simualtor, see INSTALL

2. To run a simulation, examples files can be found in the directories within the directory: sample_inputs
- to  initiate a simulations in terminal, perform "./nerdss -f parms.inp" which will call 'nerdss' to initialize '-f' a simulation by reading setup file 'parms.inp'. The simulation will create particles positional trajectory, historgram of aggregation, and simulation state file over time. 

3. To continue the simualtion,  perform "./nerdss -r restart.dat" which will resum the simulation. 


