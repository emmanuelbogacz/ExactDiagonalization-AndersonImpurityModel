# Exact diagonalization - Anderson impurity model
This repository outlines the calculation of the spectral function of the Anderson impurity model using exact diagonalization, for any temperature. 
This code is intended for educational purposes, and is presented in two flavors:
- Diagonalization of the Hamiltonian matrix (no quantum number block sectors)
- Block diagonalization of the Hamiltonian matrix according to quantum numbers $(Q,S_z)$.
The first one is more straightforward to understand and implement, whereas the latter one is used in practice since it yields better performance. Details on the Hamiltonian parameters, basis convention, as well as impurity occupation and spectral function can be found in this [pdf file](ExactDiagonalization-AndersonImpurityModel.pdf).
This code utilizes the "unfolded" geometry of the fermionic chain (https://doi.org/10.1103/PhysRevB.80.165117), which for exact diagonalization ensures that no negative signs have to be taken into account while computing the Hamiltonian tunneling term.
  
## Installation
This code is written in Julia, as it strikes a good balance between readability and performance. Installation of the julia language can be done by following the instructions on https://julialang.org/install/. This code only has two Julia dependencies, **LinearAlgebrea.jl** (part of the base package) and **Plots.jl**, which can be installed via the Julia command `import Pkg; Pkg.add("Plots")`. This package is optional and used for plotting purposes only.

## Usage
Clone this directory and edit the model parameters at the bottom of an **ed_main_\*.jl** file. This file calls the function file **ed_\*functions.jl** for each diagonalization type. The output picture **spectral_function.png** contains a plot of the spectral function and core model parameters in the title. The output data file **spectral_function.dat** is also saved for manual plotting, which contains the frequency grid as first column and the spectral function as second column.
