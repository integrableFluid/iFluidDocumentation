# iFluidCore

This class implements the general TBA equations. Following the hydrodynamical principle, the system is always in a quasi-stationary state, whereby all the methods of the class can be applied at any time for any given filling function.

Note, the `iFluidCore` class is abstract and must be extended by classes encoding an actual integrable model. [Click here to see how to implement your own model.](model.md) 


## Constructor

### `obj = iFluidCore(x_grid, rapid_grid, rapid_w, couplings, Ntypes, Options)`  
Construct an iFluidCore object containing all the information and methods of a TBA of a given model.  
**Inputs:**

- `x_grid`: Vector of `M` gridpoints in position space.  
- `rapid_grid`: Vector of `N` gridpoints in rapidity space.  
- `rapid_w`: Vector of `N` rapidity quadrature weights.  
- `couplings`: Cell array of couplings and their time and space derivatives.  
- `Ntypes`: Number of quasiparticle types in the TBA of the model.   
- `Options`: Struct of settings.   

**Returns:**

 - `obj`: `iFluidCore` object.

## Abstract (model specific) methods

### `ebare = getBareEnergy(obj, t, x, rapid, type)`
Calculate the one-particle energy.  
**Inputs:**

- `t`: Scalar indicating the time.
- `x`: Scalar or vector indicating the position.  
- `rapid`: Scalar or vector indicating the rapidity. 
- `type`: Scalar or vector indicating the type index of the quasiparticles. 

**Returns:**

 - `ebare`: One-particle energy.

---

### `pbare = getBareMomentum(obj, t, x, rapid, type)`
Calculate the one-particle momentum.  
**Inputs:**

- `t`: Scalar indicating the time.
- `x`: Scalar or vector indicating the position.  
- `rapid`: Scalar or vector indicating the rapidity. 
- `type`: Scalar or vector indicating the type index of the quasiparticles. 

**Returns:**

 - `pbare`: One-particle momentum.

---

### `de = getEnergyRapidDeriv(obj, t, x, rapid, type)`
Calculate the derivative of the one-particle energy w.r.t. the rapidity.  
**Inputs:**

- `t`: Scalar indicating the time.
- `x`: Scalar or vector indicating the position.  
- `rapid`: Scalar or vector indicating the rapidity. 
- `type`: Scalar or vector indicating the type index of the quasiparticles. 

**Returns:**

 - `de`: Derivative of one-particle energy.

---

### `dp = getMomentumRapidDeriv(obj, t, x, rapid, type)`
Calculate the derivative of the one-particle momentum w.r.t. the rapidity.  
**Inputs:**

- `t`: Scalar indicating the time.
- `x`: Scalar or vector indicating the position.  
- `rapid`: Scalar or vector indicating the rapidity. 
- `type`: Scalar or vector indicating the type index of the quasiparticles. 

**Returns:**

 - `dp`: Derivative of one-particle momentum.

---

### `dT = getScatteringRapidDeriv(obj, t, x, rapid1, rapid2, type1, type2)`
Calculate the derivative of the two-body scattering phase w.r.t. the rapidity.  
**Inputs:**

- `t`: Scalar indicating the time.
- `x`: Scalar or vector indicating the position.  
- `rapid1`: Scalar or vector indicating the (main) rapidity.
- `rapid2`: Scalar or vector indicating the (convilution) rapidity. 
- `type1`: Scalar or vector indicating the (main) type index of the quasiparticles. 
- `type2`: Scalar or vector indicating the (convolution) type index of the quasiparticles. 

**Returns:**

 - `dT`: Derivative of the two-body scattering phase. Must be an `iFluidTensor`. 

---

### `de = getEnergyCouplingDeriv(obj, coupIdx, t, x, rapid, type)`
Calculate the derivative of the one-particle energy w.r.t. the couplings.  
**Inputs:**

- `coupIdx`: Index indicating the coupling in question.
- `t`: Scalar indicating the time.
- `x`: Scalar or vector indicating the position.  
- `rapid`: Scalar or vector indicating the rapidity. 
- `type`: Scalar or vector indicating the type index of the quasiparticles. 

**Returns:**

 - `de`: Derivative of one-particle energy.

---

### `dp = getMomentumCouplingDeriv(obj, coupIdx, t, x, rapid, type)`
Calculate the derivative of the one-particle momentum w.r.t. the couplings.  
**Inputs:**

- `coupIdx`: Index indicating the coupling in question.
- `t`: Scalar indicating the time.
- `x`: Scalar or vector indicating the position.  
- `rapid`: Scalar or vector indicating the rapidity. 
- `type`: Scalar or vector indicating the type index of the quasiparticles. 

**Returns:**

 - `dp`: Derivative of one-particle momentum.

---

### `dT = getScatteringCouplingDeriv(obj, coupIdx, t, x, rapid1, rapid2, type1, type2)`
Calculate the derivative of the two-body scattering phase w.r.t. the couplings.  
**Inputs:**

- `coupIdx`: Index indicating the coupling in question.
- `t`: Scalar indicating the time.
- `x`: Scalar or vector indicating the position.  
- `rapid1`: Scalar or vector indicating the (main) rapidity.
- `rapid2`: Scalar or vector indicating the (convilution) rapidity. 
- `type1`: Scalar or vector indicating the (main) type index of the quasiparticles. 
- `type2`: Scalar or vector indicating the (convolution) type index of the quasiparticles. 

**Returns:**

 - `dT`: Derivative of the two-body scattering phase. Must be an `iFluidTensor`. 

---

## Accessor methods

### `setCouplings(obj, couplings)`  
Set the couplings of the model.  
**Inputs:**

 - `couplings`: Cell array of couplings and their time and space derivatives.  


---

### `couplings = getCouplings(obj)`  
Get the couplings of the model.  
**Returns:**

 - `couplings`: Cell array of couplings and their time and space derivatives.  


---

### `[x_grid, rapid_grid, type_grid, rapid_w] = getGrids(obj)`
Get all grids of the system.  
**Returns:**

 - `x_grid`: Vector of `M` gridpoints in position space.  
 - `rapid_grid`: Vector of `N` gridpoints in rapidity space.  
 - `rapid_w`: Vector of `N` rapidity quadrature weights.  
 - `type_grid`: Vector from 1 to `Ntypes`.  


## TBA methods


### `Q_dr = applyDressing(obj, Q, theta, t)`
Dress the quantity `Q`, threby taking into account the collective iteractions of the quasiparticles.  
**Inputs:**

- `Q`: Quantity to be dressed. Must be an `iFluidTensor`.  
- `theta`: Filling function at time `t`. Must be an `iFluidTensor`.
- `t`: Scalar indicating the time, corresponding to `theta`.

**Returns:**

 - `Q_dr`: Dressed quantity as `iFluidTensor`.

---

### `h_i = getOneParticleEV(obj, charIdx, t, x, rapid)`
Returns the one-particle eigenvalue of the i'th conserved charge. By default the 0th charge is the number operator, the 1st charge is the momentum, and the 2nd charge is the Hamiltonian.  
**Inputs:**

- `charIdx`: Scalar indicating which charge is considered.  
- `t`: Scalar indicating the time.  
- `x`: Scalar or vector indicating the position.  
- `rapid`: Scalar or vector indicating the rapidity.  

**Returns:**

 - `h_i`: Matrix containing the eigenvalues.

---

### `[rho_t, rhoS_t] = transform2rho(obj, theta_t, t_array)`
Constructs the root density and the density of states from the filling function.  
**Inputs:**

- `theta_t`: Cell array of (or single) filling functions as `iFluidTensor`.    
- `t_array`: Vector of times corresponding to the fillings in `theta_t`.   

**Returns:**

 - `rho_t`: Cell array of (or single) root densities as `iFluidTensor`.
 - `rhoS_t`: Cell array of (or single) density of states as `iFluidTensor`.

---

### `[theta_t, rhoS_t] = transform2theta(obj, rho_t, t_array)`
Constructs the filling function and the density of states from the root density.  
**Inputs:**

 - `rho_t`: Cell array of (or single) root densities as `iFluidTensor`.
 - `t_array`: Vector of times corresponding to the root densities in `rho_t`.   

**Returns:**

 - `theta_t`: Cell array of (or single) filling functions as `iFluidTensor`.
 - `rhoS_t`: Cell array of (or single) density of states as `iFluidTensor`.
---

### `[q, j] = calcCharges(obj, c_idx, theta_t, t_array)`
Calculate the expectation values of the i'th charge densities and associated currents.   
**Inputs:**

- `c_idx`: Scalar or vector of indices indicating which charges to consider.
- `theta_t`: Cell array of (or single) filling functions as `iFluidTensor`.    
- `t_array`: Vector of times corresponding to the fillings in `theta_t`.    

**Returns:**

 - `q`: Matrix of charge density expectation values for times in `t_array`.
 - `j`: Matrix of charge current expectation values for times in `t_array`.

---

### `[theta, e_eff] = calcThermalState(obj, T, TBA_couplings)`
Calculates a thermal state of the model specified by the couplings and temperature.  
**Inputs:**

- `T`: Temperature. Can be either a scalar (homogeneous) or an anonymous function of `x` (inhomogeneous).
- `TBA_couplings`: (Optional) Cell array of couplings. If none specified, the couplings of the `iFluidCore` object are used instead.  

**Returns:**

 - `theta`: Filling function of thermal state as `iFluidTensor`.
 - `e_eff`: Pseudo-energy of thermal state as `iFluidTensor`.

---

### `[v_eff, a_eff] = calcEffectiveVelocities(obj, theta, t, x, rapid, type)`
Calculate the effective velocity and acceleration of the quasiparticles given the current state of the system.  
**Inputs:**

- `theta`: Filling function at time `t`.
- `t`: Scalar indicating the time.
- `x`: Scalar or vector indicating the position.  
- `rapid`: Scalar or vector indicating the rapidity. 
- `type`: Scalar or vector indicating the type index of the quasiparticles. 

**Returns:**

 - `v_eff`: Effective velocity as `iFluidTensor`.
 - `a_eff`: Effective acceleration as `iFluidTensor`.

---

### `e_eff = calcEffectiveEnergy(obj, T, t, x, rapid)`
Calculate the pseudo-energy of the system.  
**Inputs:**

- `T`: Temperature. Can be either a scalar (homogeneous) or an anonymous function of `x` (inhomogeneous).
- `t`: Scalar indicating the time.
- `x`: Scalar or vector indicating the position.  
- `rapid`: Scalar or vector indicating the rapidity. 

**Returns:**

 - `e_eff`: Peuso-energy as `iFluidTensor`.


## Methods dependent on quasiparticle statisitcs

### `f = getStatFactor(obj, theta)`
Returns the statistical factor based the statistics followed by the quasiparticles of the model.   
**Inputs:**

 - `theta`: Filling function as `iFluidTensor`.  

**Returns:**

 - `f`: Statistical factor as `iFluidTensor`.


---

### `F = getFreeEnergy(obj, e_eff)`
Returns the free energy function of the model based on the quasiparticle statistics.  
**Inputs:**

 - `e_eff`: Pseudo-energy of the state as `iFluidTensor`.  

**Returns:**

 - `F`: Free energy function as `iFluidTensor`.


---

### `theta = calcFillingFraction(obj, e_eff)`
Calculates the filling function of a thermal state based on the quasiparticle statistics.  
**Inputs:**

 - `e_eff`: Pseudo-energy of the state as `iFluidTensor`.  

**Returns:**

 - `theta`: Filling function as `iFluidTensor`.



## Options

The `iFluidCore` class takes an `Options` struct as argument in its constructor. The `Options` struct can hold the following (case sensitive!) parameters, which are transferred to the `iFluidCore` object upon construction.  

The possible options are:

 - `homoEvol` (default `false`): Indicates if all couplings are homogeneous. If true, `a_eff` will not be calculated.
 - `tolerance` (default `1e-6`): Tolerance for iterative solution for pseudo-energy.
 - `maxcount` (default `100`): Max nnumber of iterations for finding pseudo-energy.
