# Implementing a new model

This example demonstrates how to implement a model in iFluid.

### The model

The model in question is not a real integrable model - it merely serves as a demonstration.

The implementation of a model essentially consists of programming its thermodynamic Bethe ansatz. For that we need the single particle energy and momentum along with the two-body scattering phase:  

<img src="https://latex.codecogs.com/svg.latex?\epsilon_j&space;=&space;\lambda^2&space;&plus;&space;j&space;B" title="\epsilon_j = \lambda^2 + j B" />

<img src="https://latex.codecogs.com/svg.latex?p_j&space;=&space;4&space;\lambda" title="p_j = 4 \lambda" />

<img src="https://latex.codecogs.com/svg.latex?\Theta_{i&space;j}&space;=&space;\left(&space;i&space;-&space;j&space;\right&space;)&space;\arctan&space;\left(&space;\frac{\lambda}{A}&space;\right&space;)" title="\Theta_{i j} = \left( i - j \right ) \arctan \left( \frac{\lambda}{A} \right )" />

where `j` denotes the quasiparticle type, while `A` and `B` are the couplings of our model. Furhtermore, let's assume the quasiparticles of our model to follow fermionic statistics.

### Defining the class

Create a new model by extending the `iFluidCore` class. 

```MATLAB
classdef myModel < iFluidCore
```

The constructor of our model class simply calls the constructor of the `iFluidCore` class.

```MATLAB
function obj = myModel(x_grid, rapid_grid, rapid_w, couplings, Options)   
   
    obj = obj@iFluidCore(x_grid, rapid_grid, rapid_w, couplings, Ntypes, Options);

end
```

### Implementing the TBA

Next we must implement the abstract properties and methods of the `iFluidCore` class. 
First we specify the quasipaticle statistics:

```MATLAB
properties (Access = protected)

    quasiSpecies = 'fermion'; 
    
end 
```

Next we must implement the energy, momentum and two-body scattering phase of our model along with their derivatives w.r.t. the rapidity: 

```MATLAB
function ebare = getBareEnergy(obj, t, x, rapid, type)
    % We take B to be the second coupling
    ebare = rapid.^2 + j.*obj.couplings{1,2}(t,x);
end

function pbare = getBareMomentum(obj, t, x, rapid, type)
    pbare = 4*rapid;
end

function de = getEnergyRapidDeriv(obj, t, x, rapid, type)
    de = 2*rapid;
end

function dp = getMomentumRapidDeriv(obj, t, x, rapid, type)
    % We should return something of same length as rapid 
    dp = repmat(4, length(rapid), 1);
end

function dT = getScatteringRapidDeriv(obj, t, x, rapid1, rapid2, type1, type2)
    % We take A to be the first coupling
    dT      = (type1-type2).*obj.couplings{1,1}(t,x)./( (rapid1-rapid2).^2 + obj.couplings{1,1}(t,x).^2 );
    
    % Scattering kernel should be returned as an iFluidTensor
    dT      = iFluidTensor(dT); 
end 
```

Note that all operations in general should be elementwise, as `x`, `lambda` and `type` are most often vectors.
Finally, we have to implement functions providing the derivatives w.r.t. the couplings `A` and `B`:

```MATLAB
function de = getEnergyCouplingDeriv(obj, coupIdx, t, x, rapid, type)
    if coupIdx == 1
        % Derivative w.r.t. A
        de = 0;
    else
        % Derivative w.r.t. B
        de = repmat( type, length(rapid), 1);
    end
end

    
function dp = getMomentumCouplingDeriv(obj, coupIdx, t, x, rapid, type)
    if coupIdx == 1
        % Derivative w.r.t. A
        dp = 0;
    else
        % Derivative w.r.t. B
        dp = 0;
    end
end


function dT = getScatteringCouplingDeriv(obj, coupIdx, t, x, rapid1, rapid2, type1, type2)    
    if coupIdx == 1
        % Derivative w.r.t. A
        dT = -(type1-type2).*(rapid1-rapid2)./( (rapid1-rapid2).^2 + obj.couplings{1,1}(t,x).^2 );
    else
        % Derivative w.r.t. B
        dT = 0;
    end
    
    dT = iFluidTensor(dT);
end 
```
And that's it! We have now implemented our own model, which is fully integrated in the iFluid framework.