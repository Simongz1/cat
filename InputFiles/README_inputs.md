these input files are used as a base to show the generic set up expected for each object.

several kernels where replaced by a material that computes an ADProperty -> ADKernel coupling to avoid overflow of kernel objects.

keep everything divided as: Material -> computeSomeProperty -> ADKernel -> residualContribution -> variables

IMPORTANT: ADHeatConductionTimeDerivative is replaced by MassLumpedTimeDerivative. This enforces the LHS of the Heat Equation to be dT/dt, so the RHS needs to be divided by the scaling factor $\rho_0 c_v$

for now, ignore LIPIT stuff.
