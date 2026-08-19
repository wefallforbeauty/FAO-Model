-The program should be written for real weather data or real weather data history. Then the timesfm forecast should be applicable weather data.
-Richards equation part should be much more complex: 
*A root uptake function should be added to program.
*The h parameter for Richards eq. should be adjusted to real, plausible soil statistics.
*Different numerical solutions should be investigated other than FVM like: FDM, FEM, Mixed Finite Element, Explicit Euler, Implicit Euler, Newton-Raphson, Picard.
*Adaptive time stamp could be added.
-Some uncertainty could be added to ET and \theta(x,y,z,t) functions.
-Both the forecast data and the solutions should be an input for another optimization process.
-The soil and plant data could be researched to add more realistic parameters and constants.
-Currently the ASCE-PM is taken reference for the calibration for various equations. This seems like it forces them to converge to ASCE-PM, this could be further investigated.
-A validation code should be added.
-Thornthwaite eq. should be adjusted for monthly mean temp values not for daily. This might be the reason why thornthwaite hangs at top of the graph.
-I should also investigate if the convergence of the different equations is because of any bias. Similarly relative bias could be calculated separately which a scatter plot could also help.
-The parameters and constants should be configured in such a way that they would show which equations is most susceptible to which parameter/constant and how.
-All equations, variables and numerical fluxes should be checked for dimensional consistency and unit conversions.
-Add automated water-mass conservation tests for the Richards solver.
