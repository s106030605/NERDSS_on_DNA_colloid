We modified and used EXEs/nerdss.cpp code that can read "positions.csv" as the particles' positions. 

Step1:
We initialed a random aggregation profile of 3005 particles with a approximate composition of W:R:B = 7:2:7 using NERDSS compiled from NERDSS_simulator_on_DNA_coated_colloid. 

Step2: 
We generate a position.csv for core-shell structure profile using initial_configuration_perfect_core_shell.ipynb as an example. 

Step3: 
We then restarted and continued the aggregation simulation using NERDSS compiled from NERDSS_initial_configuration. This NERDSS would read "positions.csv" to reconfigure particles' positions.
In order to create and maintain tight core-shell structure, the reaction on rate constants between particles A & B and A & C were specified as 1,000 nm^2/us in kinetics_t_stage1.csv. (The core-shell are bound together.)

Step4: 
We then restarted and continued the simulation. In order to dissolve the core and leave the ring behind, we keep on rate constants between A & C while we make the rest of off rate as 1,000s-1 in kinetics_t_stage2.csv.

 
