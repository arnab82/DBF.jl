# DBF.jl

Double Bracket Flow (DBF) algorithms for quantum Hamiltonian transformation.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://nmayhall.github.io/DBF.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://nmayhall.github.io/DBF.jl/dev/)
[![Build Status](https://github.com/nmayhall/DBF.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/nmayhall/DBF.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/nmayhall/DBF.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/nmayhall/DBF.jl)

## Overview

DBF.jl implements continuous unitary transformations of quantum Hamiltonians using the double bracket flow framework. The package provides efficient methods for:

- **Hamiltonian Diagonalization** (`dbf_diag`): Transform Hamiltonians into diagonal form
- **Ground State Energy** (`dbf_groundstate`): Find ground state energies and prepare diagonal Hamiltonians
- **Hamiltonian Disentanglement** (`dbf_disentangle`): Reduce entanglement between qubit subsystems  
- **ADAPT-VQE** (`adapt`): Adaptive variational quantum eigensolver implementation

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/nmayhall/DBF.jl")
```

## Quick Start

```julia
using DBF
using PauliOperators

# Create a 6-qubit Heisenberg Hamiltonian
N = 6
H = heisenberg_1D(N, -1.0, -1.0, -1.0)  # Jx, Jy, Jz couplings

# Prepare reference state |000000⟩
ψ = Ket{N}(0)

# Run DBF ground state algorithm
result = dbf_groundstate(H, ψ, 
    max_iter=20, 
    conv_thresh=1e-3,
    verbose=1)

# Extract final energy
E_final = expectation_value(result["hamiltonian"], ψ)
println("Ground state energy: ", E_final)
```

## Documentation

For a comprehensive guide to the theory, algorithms, and workflow, see [Theory and Workflow Documentation](docs/src/theory_and_workflow.md).

## Key Features

- **Efficient Pauli Representation**: Uses sparse Pauli operator algebra for scalability
- **Adaptive Truncation**: Dynamically removes small coefficients to control operator growth
- **Error Tracking**: Monitors truncation errors via PT2 corrections and variance
- **Flexible Generators**: Multiple generator selection strategies for different objectives
- **Built-in Hamiltonians**: Heisenberg, Hubbard, and custom models

## Methods

### Diagonalization
```julia
H_diag, generators, angles = dbf_diag(H, max_iter=100, conv_thresh=1e-7)
```

### Ground State Energy
```julia
result = dbf_groundstate(H, ψ, max_iter=20, conv_thresh=1e-3)
```

### Disentanglement
```julia
H_sep, generators, angles = dbf_disentangle(H, split_index, max_iter=100)
```

### ADAPT-VQE
```julia
pool = qubitexcitationpool(N)  # Generate operator pool
H_final, generators, angles = adapt(H, pool, ψ, max_iter=50)
```

## References

- Double Bracket Flow methodology for quantum Hamiltonian transformation
- ADAPT-VQE: [arXiv:1812.11173](https://arxiv.org/abs/1812.11173)
- Qubit Excitation Operators: [Nature Comm. Phys. 4, 1 (2021)](https://www.nature.com/articles/s42005-021-00730-0)

## License

See LICENSE file for details.
