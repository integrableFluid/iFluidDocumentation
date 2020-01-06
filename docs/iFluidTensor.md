# iFluidTensor
Quantities like the root densities of TBAs are functions of position and rapidity. Additionally, some TBAs contain multiple root densitites, each describing a different type of quasiparticle. Thus, the numerical representation of the root density can be achieved with a rank-3 tensor.  

Meanwhile, quantites like the two-body scattering phase encodes the interaction between two quasiparticles of (potentially) different type and rapidity. Thus, an extra type and rapidity arguement are needed resulting in a rank-5 tensor representation.  

Consider one of the most central operations of GHD calculations, namely the dressing of some quantity, *h*. This implies solving the set of linear equations:

<img src="https://latex.codecogs.com/svg.latex?h_{i}^{k}&space;=&space;\sum_{j}&space;\sum_{l}&space;\left(&space;\delta_{ij}&space;\delta^{kl}&space;&plus;&space;w_{j}^{l}&space;T_{ij}^{kl}&space;\vartheta_{j}^{l}&space;\right)&space;(h&space;^&space;{&space;\mathrm&space;{&space;dr&space;}&space;})_{j}^{l}&space;=&space;\sum_{j}&space;\sum_{l}&space;U_{ij}^{kl}&space;(h&space;^&space;{&space;\mathrm&space;{&space;dr&space;}&space;})_{j}^{l}" title="h_{i}^{k} = \sum_{j} \sum_{l} \left( \delta_{ij} \delta^{kl} + w_{j}^{l} T_{ij}^{kl} \vartheta_{j}^{l} \right) (h ^ { \mathrm { dr } })_{j}^{l} = \sum_{j} \sum_{l} U_{ij}^{kl} (h ^ { \mathrm { dr } })_{j}^{l}" />

We won't be going into specifics regarding the quantities of the equation, however, the subscripts mark the rapidity indices, while the superscripts mark the quasiparticle type indices. The calculation requires contracting over two indices - similar to a matrix multiplication, which is just the contraction over a single index,

The `iFluidTensor` class wraps around a standard Matlab 5-dimensional matrix and generalizes many of the standard 2-dimensional methods (such as matrix multiplication) to 4D. Thus, the multiplication above can be computed in iFluid via `h = U*h_dr` or solving the equation via `h_dr = U\h`.  

To achieve this, the `iFluidTensor` class enforces a strict convention for indices, namely  

  1. Main rapidity index.  
  2. Position index.  
  3. Main type index.  
  4. Secondary rapidity index.  
  5. Secondary type index.  

The purpose of the `iFluidTensor` is to perform multi-index contractions using the standard Matlab notation, which is achieved by overloading the mathematical operators `+  -  ./  *  .*  \  .^`.


## Constructor
An `iFluidTensor` object can be initialized in two different ways: Either via a matrix or an initializer list.

### `obj = iFluidTensor(matrix)`  
Construct an `iFluidTensor` object using the input matrix as datastructure.  
**Inputs:**

- `matrix`: Standard matlab matrix up to 5 dimensions.  

**Returns:**

 - `obj`: `iFluidTensor` object.

---

### `obj = iFluidTensor(d1, d2, ... , dN, filling )`  
Construct an `iFluidTensor` object with indices of specified sizes (up to N = 5 indices).  
**Inputs:**

- `d1, d2, ... , dN`: Integers indicating size of the *n*'th index.
- `filling`: (optional) String indicating numerical value of entries. Possible arguments are:
    - `'zeros'`: (default) All entries are 0.
    - `'ones'`: All entries are 1.
    - `'eye'`: Entries where both the two rapidity indices and the two type indices are equal are 1. All other entries are 0. (Note, tensor must be 5-dimensional).

**Returns:**

 - `obj`: `iFluidTensor` object.


## Accessor methods

### `out = getX(obj, idx, t_str)`  
Returns tensor at given position index. Equavalent to A(:,idx,:,:,:) for a standard Matlab matrix.    
**Inputs:**

- `idx`: Integer, or vector of integers, indicating the index in quation.
- `t_str`: (optional) String indicating output type:  
    - `'d'`/`'double'`: Returns standard Matlab matrix of doubles.
    - `'t'`/`'tensor'`: Returns `iFluidTensor`.

**Returns:**

- `out`: Tensor at index `idx`.

---

### `out = getRapid(obj, idx, t_str)`  
Returns tensor at given (main) rapidity index. Equavalent to A(idx,:,:,:,:) for a standard Matlab matrix.  
**Inputs:**

- `idx`: Integer, or vector of integers, indicating the index in quation.
- `t_str`: (optional) String indicating output type:  
    - `'d'`/`'double'`: Returns standard Matlab matrix of doubles.
    - `'t'`/`'tensor'`: Returns `iFluidTensor`.

**Returns:**

- `out`: Tensor at index `idx`.

---

### `out = getType(obj, idx, t_str)`  
Returns tensor at given (main) type index. Equavalent to A(:,:,idx,:,:) for a standard Matlab matrix.  
**Inputs:**

- `idx`: Integer, or vector of integers, indicating the index in quation.
- `t_str`: (optional) String indicating output type:  
    - `'d'`/`'double'`: Returns standard Matlab matrix of doubles.
    - `'t'`/`'tensor'`: Returns `iFluidTensor`.

**Returns:**

- `out`: Tensor at index `idx`.

---

### `out = getSecRapid(obj, idx, t_str)`  
Returns tensor at given (secondary) rapidity index. Equavalent to A(:,:,:,idx,:) for a standard Matlab matrix.  
**Inputs:**

- `idx`: Integer, or vector of integers, indicating the index in quation.
- `t_str`: (optional) String indicating output type:  
    - `'d'`/`'double'`: Returns standard Matlab matrix of doubles.
    - `'t'`/`'tensor'`: Returns `iFluidTensor`.

**Returns:**

- `out`: Tensor at index `idx`.

---

### `out = getSecType(obj, idx, t_str)`  
Returns tensor at given (secondary) type index. Equavalent to A(:,:,:,:,idx) for a standard Matlab matrix.  
**Inputs:**

- `idx`: Integer, or vector of integers, indicating the index in quation.
- `t_str`: (optional) String indicating output type:  
    - `'d'`/`'double'`: Returns standard Matlab matrix of doubles.
    - `'t'`/`'tensor'`: Returns `iFluidTensor`.

**Returns:**

- `out`: Tensor at index `idx`.

---

### `S = size(obj, idx)`  
Returns size of entire tensor or one of its indices.  
**Inputs:**

- `idx`: (optional) Index in question.

**Returns:**

- `S`: Sizes of all or specified index.

---

### `out = double(obj)`  
Returns tensor as standard Matlab matrix of doubles.  

**Returns:**

- `out`: Underlying matrix.



## Index manipulation methods

### `out = flatten(obj)`  
Returns a 3-dimensional matrix. The first index is the main rapidity index merged with the main type index. The second index is the secondary rapidity index merged with the secondary type index. The third index is the spatial index.  

**Returns:**

- `out`: Flattened tensor.

---

### `out = unflatten(obj)`  
Undoes the `flatten()` transformation.

**Returns:**

- `out`: Un-flattened tensor.

---

### `out = transpose(obj)` / `out = t(obj)`  
Returns tensor with both the two rapidity and the two type indices permuted.  

**Returns:**

- `out`: Transposed tensor.


## Mathematical operators
The `iFluidTensor` class overloads the elementwise mathematically operations `+  -  ./  .*  .^` along with  

  - `abs()`  
  - `log()`  
  - `exp()`  

Below are the overloaded methods, which have a different function than their Matlab counterpart.

---

### `C = mtimes(A, B)` / `C = A*B`  
Performs the two-index contraction  
<img src="https://latex.codecogs.com/svg.latex?C_{ik}^{ln}&space;=&space;\sum_{j}&space;\sum_{m}&space;A_{ij}^{lm}&space;B_{jk}^{mn}" title="C_{ik}^{ln} = \sum_{j} \sum_{m} A_{ij}^{lm} B_{jk}^{mn}" />

**Inputs:**

- `A`: Matrix or tensor with dimensions *( Ra, Xa, Ta, N, L )*, with position index length, *Xa*, equal to 1 or *M*.
- `B`: Matrix or tensor with dimensions *( N, Xb, L, Rb, Tb )*, with position index length, *Xb*, equal to 1 or *M*.

**Returns:**

- `C`: Matrix or tensor with dimensions *( Ra, max(Xa, Xb), Ta, Rb, Tb )*.

---


### `x = mldivide(A, B)` / `x = A\B`  
Solves the system of linear equations  
<img src="https://latex.codecogs.com/svg.latex?\sum_{j}&space;\sum_{m}&space;A_{ij}^{lm}&space;X_{jk}^{mn}&space;=&space;B_{ik}^{ln}" title="\sum_{j} \sum_{m} A_{ij}^{lm} X_{jk}^{mn} = B_{ik}^{ln}" />

**Inputs:**

- `A`: Matrix or tensor with dimensions *( N, Xa, L, Ra, Ta )*, with position index length, *Xa*, equal to 1 or *M*.
- `B`: Matrix or tensor with dimensions *( N, Xb, L, Rb, Tb )*, with position index length, *Xb*, equal to 1 or *M*.

**Returns:**

- `x`: Matrix or tensor with dimensions *( Ra, max(Xa, Xb), Ta, Rb, Tb )*.

---


### `B = inv(A)`  
Returns the solution to the equations  
<img src="https://latex.codecogs.com/svg.latex?\left(&space;A&space;\:&space;B&space;\right&space;)_{ij}^{lm}&space;=&space;\delta_{ij}&space;\delta^{lm}" title="\left( A \: B \right )_{ij}^{lm} = \delta_{ij} \delta^{lm}" />  

**Inputs:**

- `A`: Tensor with both rapidity indices and both type indices of equal length.

**Returns:**

- `B`: "Inverse" of `A`.

---

### `out = sum(obj, idx, t_str)`  
Returns the tensor summed over the given index.  
**Inputs:**

- `idx`: Integer, or vector of integers, indicating the index in quation.
- `t_str`: (optional) String indicating output type:  
    - `'d'`/`'double'`: Returns standard Matlab matrix of doubles.
    - `'t'`/`'tensor'`: Returns `iFluidTensor`.

**Returns:**

- `out`: Sum over `index`.
