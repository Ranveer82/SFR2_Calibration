The python modules can be used to calibrate the SFR2 parameters of MODFLOW Models develped with GMS softwatre.

Optimization Algorithm: PSO (pyswarms)

Preparation:

1. Add gages in the model at the stream observation locations.
2. create gage obrvation files as given sample and place it in the MODFLOW folder.
3. copy the MODLFOW executable into the MODFLOW foloder.
4. Also create a data file for the output time series. this will be used to map the intermittent observations to calculate the objective function.

Run the optimization for desired number of iterations and number of particles.
