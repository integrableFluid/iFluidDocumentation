## Version 1.1.0

* Added class `iFluidCorrelator` for calculating dynamic Euler-scale two-point correlations.
* Added class `LinearDiffusionSolver` for solving the GHD propagation equation with linearized diffusion.
* Added class `CollisionSolver` for solving the GHD propagation equation with added collision integral. Only works for Lieb-Liniger model in dimensional crossover.
* Changed function `calcEffectiveEnergy()` to take source term as input rather than temperature.
* Arguments `x`, `t` and `rapid` is now for many functions in `iFluidCore` optional. If no variable is passed to function, the grids stored in the class will be used instead. 