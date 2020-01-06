# iFluidSolver

The purpose of the `iFluidSolver` class is to propagate the filling function in time following the hydrodynamic equation

<img src="https://latex.codecogs.com/svg.latex?\partial&space;_&space;{&space;t&space;}&space;\vartheta&space;(\lambda)&space;&plus;&space;v&space;^&space;{&space;\mathrm&space;{&space;eff&space;}&space;}&space;[\vartheta&space;(\lambda)]&space;\:&space;\partial&space;_&space;{&space;x&space;}&space;\vartheta&space;(\lambda)&space;&plus;&space;a^&space;{&space;\mathrm&space;{&space;eff&space;}&space;}[\vartheta&space;(\lambda)]&space;\:&space;\partial&space;_&space;{&space;\lambda&space;}&space;\vartheta&space;(\lambda)&space;=&space;0" title="\partial _ { t } \vartheta (\lambda) + v ^ { \mathrm { eff } } [\vartheta (\lambda)] \: \partial _ { x } \vartheta (\lambda) + a^ { \mathrm { eff } }[\vartheta (\lambda)] \: \partial _ { \lambda } \vartheta (\lambda) = 0" />

This is achieved by approximating the propagation for a single time-step `dt`.

Note, the `iFluidSolver` class is abstract and must be extended by classes encoding algorithms for performing the time-step. [Click here to see how to implement your own algorithm.](solver.md) 


## Constructor

### `obj = iFluidSolver(coreObj, Options)`  
Construct an iFluidCore object containing all the information and methods of a TBA of a given model.  
**Inputs:**

- `coreObj`: `iFluidCore` object specifying the model and problem at hand.  
- `Options`: Struct of settings.   

**Returns:**

 - `obj`: `iFluidSolver` object.


## Abstract (algorithm specific) methods

### `[theta, u, w] = initialize(obj, theta_init, u_init, w_init, t_array )`  
Calculates and stores all the required quantities used by the algorithm implemented in `step()`.  
**Inputs:**

 - `theta_init`: Initial filling function (t = 0), as `iFluidTensor`.
 - `u_init`: Initial position characteristic (t = 0), as `iFluidTensor`.
 - `w_init`: Initial rapidity characteristic (t = 0), as `iFluidTensor`.
 - `t_array`: Vector of timesteps.  

**Returns:**

- `theta`: Filling passed to first call of `step()`. Typically equal to `theta_init`.
- `u`: Position characteristic passed to first call of `step()`. Typically equal to `u_init`.
- `w`: Rapidity characteristic passed to first call of `step()`. Typically equal to `w_init`.

---

### `[theta_next, u_next, w_next] = step(obj, theta_prev, u_prev, w_prev, t, dt)`  
Propagates the filling function one timestep `dt`.
**Inputs:**

 - `theta_prev`: Filling function at time `t`, as `iFluidTensor`.
 - `u_prev`: Position characteristic at time `t`, as `iFluidTensor`.
 - `w_prev`: Rapidity characteristic at time `t`, as `iFluidTensor`.
 - `t`: Time at the start of the step.  
 - `dt`: Length of step.  

**Returns:**

 - `theta_next`: Filling function at time `t + dt`, as `iFluidTensor`.
 - `u_next`: Position characteristic at time `t + dt`, as `iFluidTensor`.
 - `w_next`: Rapidity characteristic at time `t + dt`, as `iFluidTensor`.


## Time propagation methods

### `[theta_t, u_t, w_t] = propagateTheta(obj, theta_init, t_array)`
Propagates the filling function according to the hydrodynamical equation using the user-implemented `step()` method.   
**Inputs:**

- `theta_init`: Initial filling function (t = 0), as `iFluidTensor`.
- `t_array`: Vector of starting times of each step.

**Returns:**

- `theta_t`: Cell array of filling functions for each time in `t_array`.
- `u_t`: Cell array of position characteristics for each time in `t_array`.
- `w_t`: Cell array of rapidity characteristics for each time in `t_array`.


## Additional methods

### `[theta_next, u_next, w_next] = performFirstOrderStep(obj, theta_prev, u_prev, w_prev, t, dt)`
Propagates the filling function one timestep `dt` using a first-order approximation. This method is helpful for calculating the required quantities for more sophisticated `step()` algorithms. 
**Inputs:**

 - `theta_prev`: Filling function at time `t`, as `iFluidTensor`.
 - `u_prev`: Position characteristic at time `t`, as `iFluidTensor`.
 - `w_prev`: Rapidity characteristic at time `t`, as `iFluidTensor`.
 - `t`: Time at the start of the step.  
 - `dt`: Length of step.  

**Returns:**

 - `theta_next`: Filling function at time `t + dt`, as `iFluidTensor`.
 - `u_next`: Position characteristic at time `t + dt`, as `iFluidTensor`.
 - `w_next`: Rapidity characteristic at time `t + dt`, as `iFluidTensor`.

---

### `tensor_int = interpPhaseSpace(obj, tensor_grid, rapid_int, x_int, extrapFlag)`
Interpolates a quantity defined on the grids `x_grid` and `rapid_grid` from the stored `iFluidCore` object.  
**Inputs:**

- `tensor_grid`: Quantity to interpolate. Must be `iFluidTensor`.  
- `rapid_int`: Matrix of rapidities to interpolate to.  
- `x_int`: Matrix of positions to interpolate to.  
- `extrapFlag`: Flag for extrapolation. If `true`, extrapolation is allowed. If `false`, all extrapolated values are 0.  

**Returns:**

 - `tensor_int`: `tensor_grid` interpolate to `rapid_int` and `x_int`.


## Options

The `iFluidSolver` class takes an `Options` struct as argument in its constructor. The `Options` struct can hold the following (case sensitive!) parameters, which are transferred to the `iFluidSolver` object upon construction.  

The possible options are:

 - `extrapFlag` (default `false`): Indicates if extrapolation is allowed.
 - `periodRapid` (default `false`): Indicates if rapidity is periodic.