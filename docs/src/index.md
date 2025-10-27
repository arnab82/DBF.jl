```@meta
CurrentModule = DBF
```

# DBF

Documentation for [DBF](https://github.com/nmayhall/DBF.jl).

## Overview

DBF.jl is a Julia package implementing Double Bracket Flow (DBF) algorithms for quantum Hamiltonian transformation. The package provides efficient methods for:

- **Hamiltonian Diagonalization**: Transform Hamiltonians into diagonal form
- **Ground State Energy Calculation**: Find ground state energies via continuous unitary transformations
- **Hamiltonian Disentanglement**: Reduce entanglement between qubit subsystems
- **ADAPT-VQE**: Adaptive variational quantum eigensolver implementation

## Documentation Pages

- [Theory and Workflow](theory_and_workflow.md) - Comprehensive guide to DBF theory, algorithms, and usage

## Quick Start

```julia
using DBF
using PauliOperators

# Create a Heisenberg Hamiltonian
N = 6  # Number of qubits
H = heisenberg_1D(N, -1.0, -1.0, -1.0)

# Prepare reference state
ψ = Ket{N}(0)

# Run DBF ground state algorithm
result = dbf_groundstate(H, ψ, max_iter=20, verbose=1)

# Extract results
E_final = expectation_value(result["hamiltonian"], ψ)
println("Final energy: ", E_final)
```

## API Reference

```@index
```

```@autodocs
Modules = [DBF]
```
