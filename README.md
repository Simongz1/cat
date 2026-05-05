MARLINS Code (Microstructure-Aware Reactive Lagrangian INtegrated Shock Code)
=====

Fork "ml" to create a new MOOSE-based application.

For more information see: [https://mooseframework.inl.gov/getting_started/new_users.html#create-an-app](https://mooseframework.inl.gov/getting_started/new_users.html#create-an-app)
Code developed and maintained by [Simon Gonzalez](mailto:gonz1075@purdue.edu), [Koslowski Group, Purdue University](https://koslowskigroup.org/)
-------------------------------

The Microstructure-Aware Reactive Lagrangian INtegrated Shock Code (MARLINS Code) is a MOOSE-based application that provides a series of tools to peform shock response analysis for generated or imported PBX microstructures using a Lagrangian approach with stabilized shock tracking, kinetics models, equation of state models for reactants and products, melting, and thermal transport.

To see the required modules to install and run the application, the use of the lates Linux MOOSE distribution is recomended, which can be installed using conda following these instructions: [MOOSE-Linux Install Using Conda](https://mooseframework.inl.gov/getting_started/installation/conda.html).

A pre-packaged python interface is provided to generated input files based on a series of prompts, which allow the user to controll:

1. The microstructure to load (RANDOM, GENERATED, LOADED).
2. The binder and HE properties.
3. The different simulation paramters such as reactive model configuration, scaling, mesh size, output frequency.
4. Some paramters, such as EOS, kinetics, plasticity, are preset, but can be adapted to special needs.

Prepackaged scripts to generate P-x, P-t, Hugoniot, Pop-plot, and calculate run to detonation distances are available under /pyscripts/, where prompts can be followed. 

Eulerian capabilities using the Finite Volume method base objects provided by MOOSE (see [Finite Volume MOOSE](https://mooseframework.inl.gov/finite_volumes/fv_design.html)) are currently under development and implementation. Input files for two-phase hydrodynamic shock compression of a low density cavity are available under /inputs/, and the finite volume objects (Functors, Kernels. AuxKernels) are available under the /fv*/ directories.

This application makes use of Automatic Differentiation (AD) to form the jacobian matrix of the system of nonlinear partial differential equations. The use of AD, design decissions, installation, and considerations, can be accessed [here](https://mooseframework.inl.gov/automatic_differentiation/). Some capabilities of the Lagrangian system integrate AD and non-AD features, while the Eulerian system is fully implemented using AD-only features.