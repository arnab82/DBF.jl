```@meta
CurrentModule = DBF
```

# DBF.jl Documentation

Welcome to the documentation for DBF.jl - a Julia package implementing Double Bracket Flow (DBF) methods for quantum Hamiltonians using efficient Pauli operator representations.

## What is DBF?

Double Bracket Flow (DBF) is a computational framework for transforming quantum Hamiltonians through unitary evolution guided by commutator-based flow equations. The key innovation is using Pauli operator propagation to efficiently evolve operators in the Heisenberg picture.

**Key Features:**
- **Ground State Finding:** Find ground state energies without explicit diagonalization
- **Hamiltonian Diagonalization:** Transform Hamiltonians to diagonal form
- **Subsystem Disentanglement:** Separate quantum systems into decoupled subsystems
- **ADAPT-VQE:** Build variational quantum circuits adaptively
- **Efficient Pauli Arithmetic:** Optimized representation and evolution of Pauli operators

## Quick Start

```julia
using DBF
using PauliOperators

# Create a Heisenberg Hamiltonian
N = 6
H = heisenberg_1D(N, 1.0, 1.0, 1.0)

# Find ground state
ψ = Ket{N}(0)  # Reference state |000000⟩
result = dbf_groundstate(H, ψ, n_body=2, max_iter=50)

# Extract results
E_final = result["energies"][end]
H_evolved = result["hamiltonian"]
```

## Documentation Structure

This documentation is organized into four main sections:

### [Theory and Mathematics](theory.md)

Comprehensive explanation of the mathematical foundations:
- Pauli operator formalism
- Double Bracket Flow equations
- Pauli propagation techniques
- Optimization methods
- Truncation and error analysis
- Perturbation theory corrections

### [DBF Variants Comparison](dbf_variants.md)

**NEW!** Detailed comparison of the four DBF methods:
- Driven DBF vs Canonical DBF vs Diagonalization DBF vs Disentanglement DBF
- Mathematical formulation differences
- Implementation logic comparison
- When to use each variant
- Computational cost analysis
- Code-level differences explained

### [Repository Structure](structure.md)

Complete reference for the codebase:
- Directory and file organization
- Function-by-function documentation
- Module dependencies
- Workflow examples
- Performance considerations

### [User Guide and Examples](guide.md)

Practical tutorials and examples:
- Getting started
- Basic Pauli operations
- Ground state finding
- Hamiltonian diagonalization
- Fermionic systems
- ADAPT-VQE workflows
- Quantum gates and state preparation
- Advanced subspace methods
- Tips and best practices

## Core Algorithms

DBF.jl implements several quantum algorithms:

1. **DBF Ground State (`dbf_groundstate`)**: Finds ground state energy by evolving the Hamiltonian in a reference state using the flow equation: `dH/dt = [H, [H, P]]`

2. **DBF Diagonalization (`dbf_diag`)**: Diagonalizes operators using: `dH/dt = [[H_d, H], H]`

3. **DBF Disentanglement (`dbf_disentangle`)**: Separates subsystems by minimizing inter-subsystem coupling

4. **ADAPT (`adapt`)**: Adaptive variational quantum eigensolver that builds circuits by selecting operators from a pool based on gradients

## Key Functions

### Evolution
- `evolve(O, G, θ)` - Evolve operator: `O' = exp(iθG/2) O exp(-iθG/2)`
- `evolve!(O, G, θ)` - In-place evolution

### Hamiltonians
- `heisenberg_1D(N, Jx, Jy, Jz)` - 1D Heisenberg model
- `heisenberg_2D(Nx, Ny, Jx, Jy, Jz)` - 2D Heisenberg model
- `hubbard_model_1D(L, t, U)` - 1D Fermi-Hubbard model
- `fermi_hubbard_2D(Lx, Ly, t, U)` - 2D Fermi-Hubbard model

### Utilities
- `coeff_clip!(H, thresh)` - Truncate small coefficients
- `weight_clip!(H, max_weight)` - Truncate high-weight terms
- `inner_product(O1, O2)` - Hilbert-Schmidt inner product
- `expectation_value(H, ψ)` - Compute `⟨ψ|H|ψ⟩`
- `variance(H, ψ)` - Compute variance
- `pt2(H, ψ)` - Second-order perturbation theory

### Quantum Gates
- `hadamard(O, q)` - Hadamard gate
- `cnot(O, c, t)` - CNOT gate
- `X_gate`, `Y_gate`, `Z_gate` - Pauli gates
- `S_gate`, `T_gate` - Phase gates

## Installation

```julia
using Pkg
Pkg.add("DBF")
```

Or install directly from the repository:

```julia
Pkg.add(url="https://github.com/arnab82/DBF.jl")
```

## Dependencies

DBF.jl relies on several Julia packages:
- **PauliOperators.jl**: Core Pauli operator arithmetic
- **LinearAlgebra**: Matrix operations
- **Optim.jl**: Optimization routines
- **KrylovKit.jl**: Iterative eigensolvers
- **JLD2.jl**: Data serialization for checkpointing

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

See the [repository](https://github.com/arnab82/DBF.jl) for more details.

## License

DBF.jl is released under the MIT License. See the LICENSE file for details.

## API Reference

```@index
```

```@autodocs
Modules = [DBF]
```
