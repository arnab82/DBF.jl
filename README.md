# DBF.jl

Double Bracket Flow (DBF) for Quantum Hamiltonians using Pauli Operator Propagation

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://nmayhall.github.io/DBF.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://nmayhall.github.io/DBF.jl/dev/)
[![Build Status](https://github.com/nmayhall/DBF.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/nmayhall/DBF.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/nmayhall/DBF.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/nmayhall/DBF.jl)

## Overview

DBF.jl implements Double Bracket Flow (DBF) methods for quantum Hamiltonians using efficient Pauli operator representations and propagation techniques. The package provides multiple algorithms for:

- **Ground State Finding**: Find ground state energies without explicit diagonalization
- **Hamiltonian Diagonalization**: Transform Hamiltonians to diagonal form
- **Subsystem Disentanglement**: Separate quantum systems into decoupled subsystems  
- **ADAPT-VQE**: Build variational quantum circuits adaptively

## Quick Start

```julia
using DBF
using PauliOperators

# Create a Heisenberg Hamiltonian
N = 6
H = heisenberg_1D(N, 1.0, 1.0, 1.0)

# Find ground state energy
ψ = Ket{N}(0)  # Reference state |000000⟩
result = dbf_groundstate(H, ψ, n_body=2, max_iter=50)

# Extract results
E_final = result["energies"][end]
println("Ground state energy: $E_final")
```

## Key Features

- **Efficient Pauli Evolution**: Fast propagation of Pauli operators under unitary transformations
- **Multiple DBF Variants**: Driven DBF, Canonical DBF, Diagonalization DBF, Disentanglement DBF
- **Built-in Hamiltonians**: Heisenberg models (1D/2D), Fermi-Hubbard models, central spin models
- **Quantum Gates**: CNOT, Hadamard, X/Y/Z gates, S/T gates
- **State Preparation**: Sequences for cluster states, Néel states, RVB states
- **Error Tracking**: Monitor truncation errors, variance, PT2 corrections
- **Checkpointing**: Save progress for long calculations

## Documentation

📚 **[Complete Documentation](https://nmayhall.github.io/DBF.jl/dev/)**

The documentation includes:

1. **[Theory and Mathematics](docs/src/theory.md)**: Mathematical foundations, Pauli propagation, DBF equations
2. **[DBF Variants Comparison](docs/src/dbf_variants.md)**: Detailed comparison of the 4 DBF methods with code-level explanations
3. **[Repository Structure](docs/src/structure.md)**: Complete function reference and code organization
4. **[User Guide and Examples](docs/src/guide.md)**: 10 practical examples and tutorials

## DBF Variants

### 1. Driven DBF (`dbf_groundstate`)
```julia
result = dbf_groundstate(H, ψ, n_body=2, max_iter=100)
```
**Purpose**: Find ground state energy by evolving `H` in reference state `ψ`  
**Flow**: `dH/dt = [H, [H, P]]` where `P ≈ |ψ⟩⟨ψ|`

### 2. Canonical DBF (`groundstate_diffeq`)
```julia
H_final, gens, angles = groundstate_diffeq(H, ψ, n_body=2, stepsize=0.01)
```
**Purpose**: Same as Driven DBF but using fixed stepsize (faster per iteration)  
**Flow**: `dH/dt = [H, [H, P]]` with gradient descent

### 3. Diagonalization DBF (`dbf_diag`)
```julia
H_diag, gens, angles = dbf_diag(H, max_iter=50)
```
**Purpose**: Diagonalize Hamiltonian to find full eigenspectrum  
**Flow**: `dH/dt = [[H_d, H], H]` where `H_d = diag(H)`

### 4. Disentanglement DBF (`dbf_disentangle`)
```julia
H_separated, gens, angles = dbf_disentangle(H, M, max_iter=50)
```
**Purpose**: Separate system into subsystems at split index `M`  
**Flow**: `dH/dt = [[H_P, H], H]` where `H_P` is P-space projection

See [DBF Variants Comparison](docs/src/dbf_variants.md) for detailed differences.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/arnab82/DBF.jl")
```

## Examples

### Ground State of 2D Heisenberg Model

```julia
using DBF, PauliOperators

# 3x3 lattice
H = heisenberg_2D(3, 3, 1.0, 1.0, 1.0, periodic=true)
ψ = Ket{9}(0)

result = dbf_groundstate(H, ψ, 
    n_body=2,
    max_iter=100,
    evolve_coeff_thresh=1e-8,
    verbose=1)

println("Energy per site: $(result["energies"][end] / 9)")
```

### Fermi-Hubbard Model

```julia
# 4-site 1D Hubbard (8 qubits)
H = hubbard_model_1D(4, t=1.0, U=2.0)
ψ = Ket{8}(0b10101010)  # Half-filling

result = dbf_groundstate(H, ψ, n_body=2, max_iter=50)
```

### Hamiltonian Diagonalization

```julia
H = heisenberg_1D(6, 1.0, 1.0, 1.0)
H_diag, gens, angles = dbf_diag(H, max_iter=50)

# Extract eigenvalues
eigenvalues = sort([c for (p,c) in diag(H_diag)])
println("Ground state: $(eigenvalues[1])")
println("Gap: $(eigenvalues[2] - eigenvalues[1])")
```

### ADAPT-VQE

```julia
# Create operator pool
pool = Vector{PauliBasis{6}}()
for i in 1:6, j in i+1:6
    push!(pool, PauliBasis(Pauli(6, Y=[i], X=[j])))
end

H = heisenberg_1D(6, 1.0, 1.0, 1.0)
ψ = Ket{6}(0)

H_final, gens, angles = adapt(H, pool, ψ, max_iter=20)
println("Circuit depth: $(length(gens))")
```

## Available Hamiltonians

**Spin Systems:**
- `heisenberg_1D`, `heisenberg_2D` - Heisenberg models
- `heisenberg_central_spin` - Central spin model
- `heisenberg_sparse` - Random sparse Heisenberg

**Fermionic Systems (Jordan-Wigner):**
- `hubbard_model_1D` - 1D Fermi-Hubbard
- `fermi_hubbard_2D` - 2D Fermi-Hubbard

**Operators:**
- `S2(N)` - Total spin squared
- `Sz(N)` - Total spin z-component

## Performance Tips

1. **Use appropriate truncation thresholds**:
   ```julia
   evolve_coeff_thresh = 1e-8  # Balance accuracy vs size
   ```

2. **Limit operator weight**:
   ```julia
   evolve_weight_thresh = N-2  # Prevent exponential growth
   ```

3. **Use checkpointing for long runs**:
   ```julia
   dbf_groundstate(H, ψ, checkfile="my_calculation")
   ```

4. **Enable threading**:
   ```bash
   export JULIA_NUM_THREADS=4
   julia script.jl
   ```

## Dependencies

- [PauliOperators.jl](https://github.com/nmayhall-pnnl/PauliOperators.jl) - Core Pauli arithmetic
- LinearAlgebra, Optim, KrylovKit, JLD2, OrderedCollections

## Citation

If you use DBF.jl in your research, please cite:

```bibtex
@software{dbf_jl,
  title = {DBF.jl: Double Bracket Flow for Quantum Hamiltonians},
  author = {Mayhall, Nick and contributors},
  url = {https://github.com/arnab82/DBF.jl},
  year = {2024}
}
```

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request

## License

MIT License - see LICENSE file for details.
