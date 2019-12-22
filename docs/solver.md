# Implementing a new solver

This example demonstrates how to implement a first-order stepper to solve the GHD equation:

<img src="https://latex.codecogs.com/svg.latex?\partial&space;_&space;{&space;t&space;}&space;\vartheta&space;(\lambda)&space;&plus;&space;v&space;^&space;{&space;\mathrm&space;{&space;eff&space;}&space;}&space;[\vartheta&space;(\lambda)]&space;\:&space;\partial&space;_&space;{&space;x&space;}&space;\vartheta&space;(\lambda)&space;&plus;&space;a^&space;{&space;\mathrm&space;{&space;eff&space;}&space;}[\vartheta&space;(\lambda)]&space;\:&space;\partial&space;_&space;{&space;\lambda&space;}&space;\vartheta&space;(\lambda)&space;=&space;0" title="\partial _ { t } \vartheta (\lambda) + v ^ { \mathrm { eff } } [\vartheta (\lambda)] \: \partial _ { x } \vartheta (\lambda) + a^ { \mathrm { eff } }[\vartheta (\lambda)] \: \partial _ { \lambda } \vartheta (\lambda) = 0" />

The new solver must extend the `iFluidSolver` class and implement the abstract methods  `step()` and `initialize()`.  
The equation above admidts the implicit solution

<img src="https://latex.codecogs.com/svg.latex?\vartheta\left(t^{\prime},&space;x,&space;\lambda\right)=\vartheta(t,&space;\:&space;\tilde{x}(t^{\prime},&space;t),&space;\:&space;\tilde{\lambda}&space;(t^{\prime},&space;t&space;))" title="\vartheta\left(t^{\prime}, x, \lambda\right)=\vartheta(t, \: \tilde{x}(t^{\prime}, t), \: \tilde{\lambda} (t^{\prime}, t ))" />

where the trajectories are given by

<img src="https://latex.codecogs.com/svg.latex?\tilde{x}\left(t^{\prime},&space;t\right)=x-\int_{t}^{t^{\prime}}&space;\mathrm{d}&space;\tau&space;\:&space;v_{\tau}^{\mathrm{eff}}(\tilde{x}(\tau,&space;t),&space;\tilde{\lambda}(\tau,&space;t))" title="\tilde{x}\left(t^{\prime}, t\right)=x-\int_{t}^{t^{\prime}} \mathrm{d} \tau \: v_{\tau}^{\mathrm{eff}}(\tilde{x}(\tau, t), \tilde{\lambda}(\tau, t))" />

and

<img src="https://latex.codecogs.com/svg.latex?\tilde{\lambda}\left(t^{\prime},&space;t\right)=\lambda-\int_{t}^{t^{\prime}}&space;\mathrm{d}&space;\tau&space;\:&space;a_{\tau}^{\mathrm{eff}}(\tilde{x}(\tau,&space;t),&space;\tilde{\lambda}(\tau,&space;t))" title="\tilde{\lambda}\left(t^{\prime}, t\right)=\lambda-\int_{t}^{t^{\prime}} \mathrm{d} \tau \: a_{\tau}^{\mathrm{eff}}(\tilde{x}(\tau, t), \tilde{\lambda}(\tau, t))" />

The equations for the trajectories are also implicit, as the effective velcity and acceleration depends on the state of the system at the given time.
The solver presented in this example provides a first-order approximation to the trajectories.


### The algorithm

The iFluid solvers must provide a solution for a single time-step 

<img src="https://latex.codecogs.com/svg.latex?\vartheta(t,&space;x,&space;\lambda)&space;\to&space;\vartheta(t&space;&plus;&space;dt,&space;x,&space;\lambda)" title="\vartheta(t, x, \lambda) \to \vartheta(t + dt, x, \lambda)" />

To first order we approximate the trajectories as

<img src="https://latex.codecogs.com/svg.latex?\tilde{x}(t&space;&plus;&space;dt,&space;t&space;)&space;=&space;x&space;-&space;dt&space;\:&space;v_{t}^{\mathrm{eff}}&space;(x&space;,&space;\lambda)&space;&plus;&space;O(dt^2)" title="\tilde{x}(t + dt, t ) = x - dt \: v_{t}^{\mathrm{eff}} (x , \lambda) + O(dt^2)" />

and

<img src="https://latex.codecogs.com/svg.latex?\tilde{\lambda}(t&space;&plus;&space;dt,&space;t&space;)&space;=&space;\lambda&space;-&space;dt&space;\:&space;a_{t}^{\mathrm{eff}}&space;(x&space;,&space;\lambda)&space;&plus;&space;O(dt^2)" title="\tilde{\lambda}(t + dt, t ) = \lambda - dt \: a_{t}^{\mathrm{eff}} (x , \lambda) + O(dt^2)" />



### Defining the class

Create a new solver by extending the `iFluidSolver` class. 

```MATLAB
classdef FirstOrderSolver < iFluidSolver
```

Quantities like the effective and velocity and acceleration are calculated in the `iFluidCore` class. Therefore, we must pass an object of said class to our solver upon instantiation and store it for future use. This is taken care of by passing the object to the superclass constructor:

```MATLAB
function obj = FirstOrderSolver(coreObj, Options)        
    obj = obj@iFluidSolver(coreObj, Options);
end
```

### Implementing the stepper 

Next, we implement the abstract method `step()`, which propagates our system a single time step.

```MATLAB
function [theta_next, u_next, w_next] = performFirstOrderStep(obj, theta_prev, u_prev, w_prev, t, dt)            
    [v_eff, a_eff]  = obj.coreObj.calcEffectiveVelocities(theta_prev, t, obj.x_grid, obj.rapid_grid, obj.type_grid); 
        
    x_tilde         = obj.x_grid - dt*v_eff;
    r_tilde         = obj.rapid_grid - dt*a_eff;
    
    theta_next      = obj.interpPhaseSpace(theta_prev, r_tilde, x_tilde, obj.extrapFlag);
    u_next          = obj.interpPhaseSpace(u_prev, r_tilde, x_tilde, true ); % always extrapolate u
    w_next          = obj.interpPhaseSpace(w_prev, r_tilde, x_tilde, true ); % always extrapolate u      
end
```

The method `interpPhaseSpace()` takes care of interpolating a quantity defined on the pre-specified grids to the points of the trajectories.


### Implementing the initializer

Finally, we have the implement the abstract method `initialize()`, which prepares any quantities needed for the first time step to proceed. The first order algorithm specified here does not require such quantities. Therefore, we don't have to do anything special in the `initialize()` method.

``` MATLAB
function [theta, u, w] = initialize(obj, theta_init, u_init, w_init, t_array)
    theta   = theta_init;
    u       = u_init;
    w       = w_init;
end
```

### Application

Once the steps above are completed, the newly implemented solver can be instantiated and used for solving the GHD equation:

```MATLAB
Solver1 = FirstOrderSolver( someModel );
theta_t = Solver1.propagateTheta( theta_init, t_array );
```

Related:
[How to implement you own model](model.md)