# Repository Structure and Function Reference

## Overview

This document provides a comprehensive guide to the DBF.jl repository structure, detailing what each file contains and what each function does.

## Directory Structure

```
DBF.jl/
├── src/               # Source code
│   ├── DBF.jl         # Main module definition
│   ├── helpers.jl     # Helper functions and utilities
│   ├── evolve.jl      # Pauli operator evolution
│   ├── groundstate.jl # Ground state finding algorithms
│   ├── diagonalization.jl  # Diagonalization via DBF
│   ├── disentangle.jl # Subsystem disentanglement
│   ├── adapt.jl       # ADAPT-style variational methods
│   ├── hamiltonians.jl # Hamiltonian generators
│   ├── schrodinger_picture.jl # State-based operations
│   └── test2.jl       # Experimental/testing code
├── test/              # Test suite
├── docs/              # Documentation
└── Project.toml       # Package dependencies
```

## Source Files

### src/DBF.jl

**Purpose:** Main module file that defines the DBF module and exports public functions.

**Key Components:**
- Type alias: `XZPauliSum{T}` - Efficient storage of Pauli operators indexed by X and Z bits
- Module imports: PauliOperators, Printf, LinearAlgebra, OrderedCollections
- Includes all source files
- Exports public API functions

**Exported Functions:**
- `evolve`, `evolve!` - Pauli operator evolution
- `inner_product` - Inner product of Pauli sums
- `offdiag` - Extract off-diagonal part
- `dbf_diag` - Diagonalization
- `dbf_groundstate` - Ground state finding
- `dbf_disentangle` - Subsystem disentanglement
- `adapt` - ADAPT-style optimization
- `coeff_clip!` - Truncate small coefficients
- `pack_x_z` - Convert to XZ-indexed format
- `project` - Project onto subspace
- Quantum gates: `hadamard`, `cnot`, `X_gate`, `Y_gate`, `Z_gate`, `S_gate`, `T_gate`

---

### src/helpers.jl

**Purpose:** Utility functions for Pauli operator manipulation, weight calculations, and data structures.

**Key Functions:**

#### Weight Calculations

```julia
majorana_weight(Pb::Union{PauliBasis{N}, Pauli{N}}) where N
```
Computes the Majorana fermion weight of a Pauli string. This is important for fermionic systems mapped to qubits via Jordan-Wigner transformation.

**Algorithm:**
- Iterates through Pauli string from right to left
- Tracks control bits for X operators
- Z operators contribute weight 2 when controlled, otherwise 0
- X operators flip control and contribute weight 1

```julia
pauli_weight(Pb::Union{PauliBasis{N}, Pauli{N}}) where N
```
Computes standard Pauli weight (number of non-identity operators).

#### Inner Products and Norms

```julia
inner_product(O1::PauliSum{N,T}, O2::PauliSum{N,T}) where {N,T}
```
Computes the Hilbert-Schmidt inner product: `⟨O1|O2⟩ = Tr(O1† O2)`.

```julia
norm(p::PauliSum{N,T}) where {N,T}
```
Computes Frobenius norm: `||O|| = √(Σ|c_i|²)`.

#### Operator Manipulation

```julia
largest_diag(ps::PauliSum{N,T}) where {N,T}
```
Finds the largest diagonal term.

```julia
diag(ps::PauliSum{N,T}) where {N,T}
```
Extracts diagonal part (terms with no X operators).

```julia
offdiag(ps::PauliSum{N,T}) where {N,T}
```
Extracts off-diagonal part (terms with X operators).

```julia
weight(p::PauliBasis)
```
Returns Pauli weight using bitwise operations: `count_ones(p.x | p.z)`.

#### Truncation Functions

```julia
coeff_clip!(ps::PauliSum{N}; thresh=1e-16) where {N}
```
Removes terms with `|coefficient| < thresh`.

```julia
weight_clip!(ps::PauliSum{N}, max_weight::Int) where {N}
```
Removes terms with Pauli weight > `max_weight`.

```julia
majorana_weight_clip!(ps::PauliSum{N}, max_weight::Int) where {N}
```
Removes terms with Majorana weight > `max_weight`.

#### Advanced Functions

```julia
find_top_k(dict, k=10)
```
Efficiently finds the k largest terms by absolute value. Optimized for `k << length(dict)`.

```julia
find_top_k_offdiag(dict, k=10)
```
Same as `find_top_k` but only considers off-diagonal terms.

```julia
get_weight_counts(O::PauliSum{N}) where N
```
Returns histogram of term counts by weight.

```julia
get_weight_probs(O::PauliSum{N}) where N
```
Returns probability distribution `P(w) = Σ_{weight(p)=w} |c_p|²`.

#### Matrix Representations

```julia
Matrix(O::PauliSum{N,T}, S::Vector{Ket{N}}) where {N,T}
```
Builds matrix representation in subspace defined by basis `S`.

**Algorithm:**
1. Converts to XZ-indexed format for efficiency
2. Computes diagonal elements
3. Computes off-diagonal elements using XOR to find connected states
4. Returns Hermitian matrix

```julia
Matrix(O::XZPauliSum{T}, basis::Vector{Ket{N}}) where {N,T}
```
Optimized version for pre-packed operators.

#### Data Structure Conversion

```julia
pack_x_z(H::PauliSum{N,T}) where {N,T}
```
Converts `PauliSum` to `XZPauliSum` format: `Dict{X_bits, Vector{(Z_bits, coeff)}}`.

**Purpose:** Enables efficient matrix-vector products by grouping terms with same X pattern.

---

### src/evolve.jl

**Purpose:** Core evolution functions for Pauli operators and quantum states under unitary rotations.

**Key Functions:**

#### Operator Evolution

```julia
evolve(O::PauliSum{N, T}, G::PauliBasis{N}, θ::Real) where {N,T}
```
Evolves operator: `O(θ) = exp(iθG/2) O exp(-iθG/2)`.

**Algorithm:**
1. Separates commuting and non-commuting terms
2. Commuting terms: multiplied by `cos(θ)`
3. Non-commuting terms: generates `sin(θ)` contribution via `G×O`
4. Returns `O cos(θ) - i sin(θ) [G,O]`

```julia
evolve!(O::PauliSum{N, T}, G::PauliBasis{N}, θ::Real) where {N,T}
```
In-place version of `evolve`.

```julia
evolve(O0::PauliSum{N,T}, g::Vector{PauliBasis{N}}, θ::Vector{<:Real}; kwargs...) where {N,T}
```
Applies sequence of rotations with truncation and error tracking.

**Parameters:**
- `thresh` - coefficient truncation threshold
- `max_weight` - maximum Pauli weight
- `verbose` - verbosity level
- `compute_var_error` - track variance errors
- `ψ` - reference state for expectation values

#### State Evolution

```julia
evolve(K::KetSum{N, T}, G::PauliBasis{N}, θ::Real) where {N,T}
```
Evolves quantum state: `|ψ(θ)⟩ = exp(-iθG/2)|ψ⟩`.

**Algorithm:**
1. Apply `cos(θ/2)` to original state
2. Apply `G` to state and multiply by `i sin(θ/2)`
3. Combine contributions

#### Dissipation

```julia
dissipate!(O::PauliSum, lmax::Int, γ::Real)
```
Applies exponential damping to high-weight terms: `c_p → exp(-γ(w-lmax)) c_p` for `w > lmax`.

#### Quantum Gates

```julia
cnot(p::Union{PauliSum{N}, KetSum{N}}, c::Int, t::Int) where N
```
Applies CNOT gate using Pauli rotations:
1. Rotate by `ZX_ct` (π/2)
2. Rotate by `X_t` (-π/2)
3. Rotate by `Z_c` (-π/2)

```julia
hadamard(p::Union{PauliSum{N}, KetSum{N}}, q::Int) where N
```
Applies Hadamard gate using rotations:
1. Z rotation (π/2)
2. X rotation (π/2)
3. Z rotation (π/2)

```julia
X_gate(p::Union{PauliSum{N}, KetSum{N}}, q) where N
Y_gate(p::Union{PauliSum{N}, KetSum{N}}, q) where N
Z_gate(p::Union{PauliSum{N}, KetSum{N}}, q) where N
S_gate(p::Union{PauliSum{N}, KetSum{N}}, q) where N
T_gate(p::Union{PauliSum{N}, KetSum{N}}, q) where N
```
Standard single-qubit gates implemented via Pauli evolution.

#### State Preparation Sequences

```julia
get_1d_neel_state_sequence(N)
```
Returns gate sequence for 1D Néel state (alternating ↑↓↑↓...).

```julia
get_1d_cluster_state_sequence(N)
```
Returns gate sequence for 1D cluster state.

```julia
get_rvb_sequence(N)
```
Returns gate sequence for RVB (Resonating Valence Bond) state.

**Sequences return:** `(generators::Vector{PauliBasis{N}}, angles::Vector{Float64})`

---

### src/groundstate.jl

**Purpose:** Ground state finding via Double Bracket Flow.

**Key Functions:**

#### Main Algorithm

```julia
dbf_groundstate(Oin::PauliSum{N,T}, ψ::Ket{N}; kwargs...) where {N,T}
```
Finds ground state energy by evolving Hamiltonian in reference state.

**Parameters:**
- `n_body` - order of projector approximation (default: 1)
- `max_iter` - maximum iterations
- `conv_thresh` - convergence threshold on gradient norm
- `evolve_coeff_thresh` - truncation threshold during evolution
- `evolve_weight_thresh` - max weight during evolution
- `evolve_mweight_thresh` - max Majorana weight during evolution
- `grad_coeff_thresh` - threshold for gradient terms
- `grad_weight_thresh` - max weight for gradient
- `energy_lowering_thresh` - minimum energy decrease to accept rotation
- `max_rots_per_grad` - max rotations per gradient computation
- `clifford_check` - check for Clifford rotations (π/2)
- `compute_var_error` - track variance errors
- `compute_pt2_error` - track PT2 errors
- `checkfile` - checkpoint file path

**Algorithm:**
1. Create reference projector `P = Σ Z_i Z_j ... Z_k` (n-body)
2. Compute gradient: `G = [P, H]`
3. Truncate gradient by coefficient and weight
4. For each significant gradient term:
   - Optimize angle `θ` to minimize `⟨ψ|H(θ)|ψ⟩`
   - Evolve `H`
   - Truncate evolved operator
   - Track energy and errors
5. Repeat until gradient norm < threshold

**Returns:**
- Dictionary with:
  - `hamiltonian` - final evolved Hamiltonian
  - `energies` - energy at each rotation
  - `variances` - variance at each rotation
  - `accumulated_error` - total truncation error
  - `generators` - list of generators used
  - `angles` - list of angles used
  - `pt2_per_grad` - PT2 correction at each iteration

#### Helper Functions

```julia
optimize_theta_expval(O::PauliSum{N,T}, G::PauliBasis{N}, ψ::Ket{N}; verbose=1) where {N,T}
```
Finds optimal rotation angle to minimize energy.

**Cost Function:**
```math
E(θ) = cos²(θ/2) ⟨ψ|H|ψ⟩ + sin²(θ/2) ⟨Gψ|H|Gψ⟩ - 2i cos(θ/2)sin(θ/2) ⟨ψ|H|Gψ⟩
```

Returns `(θ_optimal, cost_function)`.

```julia
create_0_projector(N, n_body)
```
Creates n-body approximation to `|00...0⟩⟨00...0|`:
```math
P ≈ Σ_{i_1 < ... < i_k, k≤n} Z_{i_1} ⊗ ... ⊗ Z_{i_k}
```

```julia
commutator(O1::PauliSum{N}, O2::PauliSum{N}; thresh=1e-12) where N
```
Computes `[O1, O2] = O1 O2 - O2 O1` with truncation.

---

### src/diagonalization.jl

**Purpose:** Diagonalization of operators via Double Bracket Flow.

**Key Functions:**

```julia
dbf_diag(Oin::PauliSum{N,T}; kwargs...) where {N,T}
```
Diagonalizes operator by maximizing diagonal norm.

**Parameters:**
- `max_iter` - maximum iterations
- `thresh` - coefficient threshold
- `conv_thresh` - convergence threshold
- `evolve_coeff_thresh` - truncation during evolution
- `evolve_weigth_thresh` - max weight during evolution
- `search_n_top` - number of top terms to consider in commutator

**Algorithm:**
1. Compute diagonal part: `S = diag(H)`
2. Find top commutator terms: `[S, H]`
3. Select generator with largest coefficient
4. Optimize angle to maximize `||diag(H(θ))||`
5. Evolve and truncate
6. Repeat until converged

**Returns:** `(H_diag, generators, angles)`

#### Optimization

```julia
optimize_theta_diagonalization(O, G; verbose=1)
```
Finds angle maximizing diagonal norm using analytical formula (see theory document).

#### Commutator Utilities

```julia
max_of_commutator2(A::PauliSum{N}, B::PauliSum{N}; n_top=1000) where N
```
Efficiently computes commutator using only top terms from each operator.

```julia
diag_GOG(G::PauliBasis{N}, O::PauliSum{N}) where N
```
Computes diagonal part of `G O G`.

```julia
diag_commutator(G::PauliBasis{N}, O::PauliSum{N}) where N
```
Computes diagonal part of `[G, O]`.

---

### src/disentangle.jl

**Purpose:** Separate quantum system into two subsystems (P-space and Q-space).

**Key Functions:**

```julia
dbf_disentangle(Oin::PauliSum{N,T}, M::Int; kwargs...) where {N,T}
```
Disentangles system at qubit index `M`, minimizing Q-space (qubits > M) coupling.

**Parameters:**
- `M` - split index between P and Q spaces
- `max_iter`, `thresh`, `conv_thresh`, `evolve_coeff_thresh`

**Algorithm:**
1. Project to P-space: `S = P(H)`
2. Find generator from `[S, H]`
3. Optimize to minimize `||Q(H(θ))||`
4. Evolve and repeat

```julia
p_space(O::PauliSum{N,T}, split_idx) where {N,T}
```
Extracts P-space operators (acting on qubits ≤ split_idx).

```julia
q_space(O::PauliSum{N,T}, split_idx) where {N,T}
```
Extracts Q-space operators (acting across split).

```julia
optimize_theta_disentangle(O, G, M; stepsize=.001, verbose=1)
```
Optimizes angle to minimize Q-space norm.

---

### src/adapt.jl

**Purpose:** ADAPT-style variational quantum eigensolver.

**Key Functions:**

```julia
adapt(Oin::PauliSum{N,T}, pool::Vector{PauliBasis{N}}, ψ::Ket{N}; kwargs...) where {N,T}
```
Builds variational circuit by selecting generators from pool.

**Parameters:**
- `pool` - available Pauli operators
- `ψ` - reference state
- `max_iter` - maximum iterations
- `conv_thresh` - convergence on gradient norm
- `grad_coeff_thresh` - threshold to include gradient term
- Additional truncation parameters

**Algorithm:**
1. For each operator in pool, compute gradient:
   ```math
   ∇_i = Im(⟨G_i ψ|H|ψ⟩ - ⟨ψ|H|G_i ψ⟩)
   ```
2. Sort by gradient magnitude
3. For significant gradients:
   - Optimize angle
   - Accept if energy lowering > threshold
   - Evolve and truncate
4. Repeat until convergence

**Returns:** `(H_final, generators, angles)`

#### Variance and Related

```julia
variance(O::PauliSum{N}, ψ::Ket{N}) where N
```
Computes `Var(H) = ⟨H²⟩ - ⟨H⟩²`.

```julia
skewness(O::PauliSum{N,T}, ψ::Ket{N}) where {N,T}
```
Computes third cumulant `κ_3 = ⟨H³⟩ - 3⟨H²⟩⟨H⟩ + 2⟨H⟩³`.

```julia
entropy(O)
```
Computes Shannon entropy: `S = -Σ_i p_i log(p_i)` where `p_i = |c_i|²/||O||²`.

#### Pool Generators

```julia
generate_pool_1_weight(N)
generate_pool_2_weight(N)
generate_pool_3_weight(N)
generate_pool_4_weight(N)
generate_pool_5_weight(N)
generate_pool_6_weight(N)
```
Generate pools of Pauli operators with specified weights.

```julia
qubitexcitationpool(n_system::Int)
```
Generates qubit excitation pool based on "Communications Physics 4, 1 (2021)".
Includes singles and doubles excitations.

---

### src/hamiltonians.jl

**Purpose:** Library of common quantum Hamiltonians.

**Spin Systems:**

```julia
heisenberg_1D(N, Jx, Jy, Jz; x=0, y=0, z=0)
```
1D Heisenberg chain with periodic boundary conditions:
```math
H = Σ_i [-2J_x X_i X_{i+1} - 2J_y Y_i Y_{i+1} - 2J_z Z_i Z_{i+1}] + Σ_i [x X_i + y Y_i + z Z_i]
```

```julia
heisenberg_2D(Nx, Ny, Jx, Jy, Jz; x=0, y=0, z=0, periodic=true)
```
2D Heisenberg model on square lattice.

```julia
heisenberg_2D_zigzag(Nx, Ny, Jx, Jy, Jz; x=0, y=0, z=0, periodic=true)
```
2D Heisenberg with zigzag (snake) indexing.

```julia
af_heisenberg(N, Jx, Jy, Jz; x=0, y=0, z=0, periodic=true)
```
Antiferromagnetic Heisenberg (factor of 1/4 on couplings).

```julia
heisenberg_central_spin(N, Jx, Jy, Jz; x=0, y=0, z=0, α=0, seed=1)
```
Central spin model: all spins couple to spin 1 with random disorder `α`.

```julia
heisenberg_sparse(N, Jx, Jy, Jz, sparsity; x=0, y=0, z=0, seed=1, α=1)
```
Random sparse Heisenberg with coupling probability `sparsity`.

**Fermionic Systems (via Jordan-Wigner):**

```julia
hubbard_model_1D(L::Int64, t::Float64, U::Float64)
```
1D Fermi-Hubbard model:
```math
H = -t Σ_{i,σ} (c†_{i,σ} c_{i+1,σ} + h.c.) + U Σ_i n_{i,↑} n_{i,↓}
```

```julia
fermi_hubbard_2D(Lx::Int, Ly::Int, t::Float64, U::Float64)
```
2D Fermi-Hubbard model.

```julia
fermi_hubbard_2D_zigzag(Lx::Int, Ly::Int, t::Float64, U::Float64)
```
2D Fermi-Hubbard with zigzag indexing.

**Helper for Jordan-Wigner:**

```julia
JWmapping(N; i::Int, j::Int)
```
Returns `c†_i c_j` in Pauli representation using Jordan-Wigner transformation.

**Utilities:**

```julia
S2(N)
```
Total spin squared operator.

```julia
Sz(N)
```
Total spin z-component.

```julia
graph_adjacency(O::PauliSum{N,T}) where {N,T}
```
Constructs adjacency matrix from Hamiltonian connectivity.

```julia
graph_laplacian(O::PauliSum{N,T}) where {N,T}
```
Constructs graph Laplacian from Hamiltonian.

---

### src/schrodinger_picture.jl

**Purpose:** Operations in the Schrödinger picture (state evolution, matrix-vector products).

**Key Functions:**

#### Matrix-Vector Products

```julia
matvec(O::XZPauliSum{T}, v::KetSum{N}) where {N,T}
```
Computes `O|ψ⟩` efficiently using XZ-indexed format.

```julia
subspace_matvec(O::XZPauliSum, v::KetSum{N,T}) where {N,T}
```
Computes `O|ψ⟩` restricted to the subspace spanned by `v`.

```julia
subspace_matvec_thread!(s::KetSum{N,T}, O::XZPauliSum, v::KetSum{N,T}) where {N,T}
```
Threaded version of subspace matrix-vector product.

#### Linear Maps

```julia
LinearMap(O::XZPauliSum{T}, basis::Vector{Ket{N}}; kwargs...) where {N,T}
```
Creates a `LinearMap` object for use with iterative eigensolvers (KrylovKit).

**Parameters:**
- `ishermitian` - operator is Hermitian (default: true)
- `issymmetric` - operator is symmetric (default: true)

#### Perturbation Theory

```julia
pt2(H::PauliSum{N,T}, ψ::Ket{N}) where {N,T}
```
Computes second-order perturbation theory correction:
```math
E^{(2)} = Σ_{x≠0} |⟨σ|H|ψ⟩|² / (E_0 - ⟨σ|H_d|σ⟩)
```

Returns `(E0, E2)`.

#### Advanced Methods

```julia
cepa(H::PauliSum, ref::Ket{N}; thresh=1e-4, verbose=4, x0=nothing, tol=1e-6) where N
```
Coupled Electron Pair Approximation - solves for correlation energy.

**Equations:**
```
P H P c + P H Q c = E P c
Q H P c + Q H Q c = E Q c
```

Solves: `(E Q - Q H Q) Q c = Q H P c`

```julia
fois_ci(Hin::PauliSum, ref::Ket{N}; thresh=1e-4, verbose=4, v0=nothing, tol=1e-6, max_iter=10, krylov_order=1) where N
```
First-Order Interacting Space CI - builds Krylov space and diagonalizes.

**Algorithm:**
1. Build basis: `{|ψ⟩, H|ψ⟩, H²|ψ⟩, ...}` up to `krylov_order`
2. Diagonalize `H` in this space using KrylovKit

#### Utility Functions

```julia
project(k::KetSum{N,T}, basis::Vector{Ket{N}}) where {N,T}
```
Projects state onto subspace defined by basis.

```julia
fill!(k::KetSum{N}, v::Vector{T}, basis::Vector{Ket{N}}) where {N,T}
```
Fills KetSum from vector using basis ordering.

---

## Test Files

### test/test_evolve.jl
Tests for Pauli evolution, gates, and state preparation.

### test/test_groundstate_dbf.jl
Tests for ground state finding algorithm.

### test/test_diag_dbf.jl
Tests for diagonalization.

### test/test_disentangle_dbf.jl
Tests for disentanglement.

### test/test_adapt.jl
Tests for ADAPT algorithm.

### test/test_helpers.jl
Tests for helper functions.

### test/test_majorana_weight.jl
Tests for Majorana weight calculations.

### test/test_schrodinger_picture.jl
Tests for state-based operations.

### test/test_theta_opt.jl
Tests for angle optimization.

---

## Workflow Summary

### Typical Ground State Workflow

```julia
using DBF
using PauliOperators

# Create Hamiltonian
N = 6
H = heisenberg_1D(N, 1.0, 1.0, 1.0)

# Reference state
ψ = Ket{N}(0)  # |000000⟩

# Run DBF ground state
result = dbf_groundstate(H, ψ, 
    n_body=2,
    max_iter=100,
    evolve_coeff_thresh=1e-8,
    grad_coeff_thresh=1e-6)

# Extract results
E_final = result["energies"][end]
H_evolved = result["hamiltonian"]
generators = result["generators"]
angles = result["angles"]
```

### Typical Diagonalization Workflow

```julia
# Create operator
H = heisenberg_2D(3, 3, 1.0, 1.0, 0.5)

# Diagonalize
H_diag, gens, θs = dbf_diag(H,
    max_iter=50,
    evolve_coeff_thresh=1e-10,
    search_n_top=100)

# Check diagonality
println("Diagonal norm: ", norm(diag(H_diag)))
println("Off-diagonal norm: ", norm(offdiag(H_diag)))
```

### Typical ADAPT Workflow

```julia
# Create Hamiltonian and pool
H = hubbard_model_1D(4, 1.0, 2.0)
pool = qubitexcitationpool(8)  # 4 sites × 2 spins

# Reference state (half-filling)
ψ = Ket{8}(0b10101010)

# Run ADAPT
H_final, gens, θs = adapt(H, pool, ψ,
    max_iter=20,
    grad_coeff_thresh=1e-6,
    evolve_coeff_thresh=1e-10)

# Final energy
E = expectation_value(H_final, ψ)
```

---

## Performance Considerations

1. **Truncation Thresholds:** Balance accuracy vs. operator size
   - `evolve_coeff_thresh`: typically 1e-8 to 1e-12
   - `grad_coeff_thresh`: typically 1e-6 to 1e-8

2. **Weight Limits:** Prevent exponential growth
   - Pauli weight: typically N/2 to N
   - Majorana weight: often more permissive

3. **Search Pool Size:** For diagonalization
   - `search_n_top`: 100-1000 depending on system size

4. **Threading:** Use `JULIA_NUM_THREADS` for parallel operations

5. **Memory:** Operator size grows; use checkpointing for large systems

---

## Dependencies

- **PauliOperators.jl:** Core Pauli arithmetic
- **LinearAlgebra:** Matrix operations
- **Optim.jl:** Optimization routines
- **OrderedCollections:** Efficient dictionaries
- **JLD2:** Checkpointing
- **TimerOutputs:** Performance profiling
- **KrylovKit:** Iterative eigensolvers
- **LinearMaps:** Linear operator interface

---

## Extending the Package

### Adding New Hamiltonians

Add to `src/hamiltonians.jl`:
```julia
function my_hamiltonian(N, params...)
    H = PauliSum(N, Float64)
    # Build H using Pauli operators
    return H
end
```

### Adding New Pools

Add to `src/adapt.jl`:
```julia
function generate_pool_custom(N)
    pool = Vector{PauliBasis{N}}([])
    # Add desired Pauli operators
    return pool
end
```

### Custom Evolution Schemes

Implement custom evolution by:
1. Computing gradients
2. Optimizing angles
3. Calling `evolve!` with generators and angles
4. Tracking errors and convergence
