# DBF Variants: Detailed Comparison

## Overview

DBF.jl implements multiple variants of the Double Bracket Flow method, each designed for different quantum computational tasks. This document explains the key differences in theory, implementation logic, and use cases.

## The Four Main DBF Variants

### 1. **Driven DBF (Ground State Finding)** - `dbf_groundstate`
### 2. **Canonical DBF (Differential Equation)** - `groundstate_diffeq`
### 3. **Diagonalization DBF** - `dbf_diag`
### 4. **Disentanglement DBF** - `dbf_disentangle`

---

## Mathematical Formulation Differences

### 1. Driven DBF (Ground State Finding)

**Flow Equation:**
```math
\frac{dH}{dt} = [H, [H, P]]
```

where `P ≈ |ψ_ref⟩⟨ψ_ref|` is the projector onto the reference state.

**Key Characteristics:**
- **Driving Operator:** Fixed reference state projector `P`
- **Goal:** Minimize `⟨ψ|H(t)|ψ⟩` (find ground state energy in reference state)
- **Optimization:** Each generator angle is optimized to minimize energy
- **Convergence:** Stops when gradient `||[P, H]||` is small

**N-body Approximation:**
```math
P ≈ \frac{1}{2^N} \sum_{i_1 < ... < i_k, k≤n} Z_{i_1} ⊗ ... ⊗ Z_{i_k}
```

This approximates the projector using up to n-body Pauli Z operators.

---

### 2. Canonical DBF (Differential Equation Approach)

**Flow Equation:**
```math
\frac{dH}{dt} = [H, [H, P]]
```

Same as Driven DBF, but with different discretization.

**Key Characteristics:**
- **Driving Operator:** Same n-body projector as Driven DBF
- **Goal:** Same as Driven DBF (minimize energy)
- **Optimization:** Uses gradient directly with fixed stepsize (no angle optimization per generator)
- **Convergence:** Stops when `||[P, H]||` is small
- **Discretization:** `θ_i = -Im(c_i) × stepsize` where `c_i` are commutator coefficients

**Key Difference from Driven DBF:**
- **Canonical:** Takes small steps along gradient direction (like gradient descent)
- **Driven:** Optimizes each rotation angle individually (slower but more accurate per rotation)

---

### 3. Diagonalization DBF

**Flow Equation:**
```math
\frac{dH}{dt} = [[H_d, H], H]
```

where `H_d = diag(H)` is the diagonal part (computational basis).

**Key Characteristics:**
- **Driving Operator:** Diagonal part `H_d` (changes at each iteration)
- **Goal:** Maximize `||diag(H(t))||` (make H as diagonal as possible)
- **Optimization:** Angle optimized to maximize diagonal norm
- **Convergence:** Stops when off-diagonal terms are negligible or generator repeats

**Cost Function:**
```math
C(θ) = ||diag(e^{iθG/2} H e^{-iθG/2})||²
```

The angle maximizes diagonal content, not minimizes energy.

---

### 4. Disentanglement DBF

**Flow Equation:**
```math
\frac{dH}{dt} = [[H_P, H], H]
```

where `H_P = P(H)` is the projection of H onto the P-space (subsystem A).

**Key Characteristics:**
- **Driving Operator:** P-space projection `P(H)` (changes at each iteration)
- **Goal:** Minimize `||Q(H(t))||` (minimize Q-space coupling between subsystems)
- **Split:** System divided at qubit index M into P-space (qubits 1 to M) and Q-space (qubits M+1 to N)
- **Optimization:** Angle optimized to minimize Q-space norm
- **Convergence:** Stops when `||Q(H)||` is small

**Space Definitions:**
- **P-space:** Operators acting only on qubits ≤ M
- **Q-space:** Operators coupling across the M boundary

---

## Comparison Table

| Aspect | Driven DBF | Canonical DBF | Diagonalization DBF | Disentanglement DBF |
|--------|-----------|---------------|-------------------|-------------------|
| **Flow** | `[H,[H,P]]` | `[H,[H,P]]` | `[[H_d,H],H]` | `[[H_P,H],H]` |
| **Driving Op** | Fixed `P` | Fixed `P` | Changing `H_d` | Changing `H_P` |
| **Objective** | Min `⟨H⟩` | Min `⟨H⟩` | Max `‖diag(H)‖` | Min `‖Q(H)‖` |
| **Gradient** | `[P, H]` | `[P, H]` | `[diag(H), H]` | `[P(H), H]` |
| **Angle Opt** | Per-generator | Fixed step | Per-generator | Per-generator |
| **Stepsize** | Optimized `θ` | `θ = c × ε` | Optimized `θ` | Optimized `θ` |
| **Reference State** | Required | Required | Not needed | Not needed |
| **Output** | Low-energy H | Low-energy H | Diagonal H | Block-diagonal H |

---

## Implementation Logic Differences

### Driven DBF Implementation

**Location:** `src/groundstate.jl:77` - function `dbf_groundstate`

**Logic Flow:**
```julia
1. Create reference projector: P = create_0_projector(N, n_body)
2. Repeat until converged:
   a. Compute gradient: G = commutator(P, O)  # [P, H]
   b. Truncate gradient by coefficient and weight
   c. Pack H into XZ format for efficient matvec
   d. Compute gradient vector: ∇_i = 2 Re⟨G_i ψ|H|ψ⟩
   e. Sort gradients by magnitude
   f. For each significant gradient:
      - Optimize angle: θ_i = argmin_θ ⟨ψ|H(θ)|ψ⟩
      - Evolve: H ← exp(iθ_i G_i/2) H exp(-iθ_i G_i/2)
      - Truncate H
      - Track energy and errors
   g. Check convergence: ||∇|| < threshold
```

**Key Code Sections:**
```julia
# Line 166: Create projector
P = create_0_projector(N, n_body)

# Line 173: Compute commutator gradient
G = commutator(P, O)

# Line 193-208: Compute gradient vector with matvec
for (p,c) in G 
    ci, σ = p*ψ
    gi = 2*real(get(σv, σ, T(0)) * c * ci)  # Gradient component
    push!(grad_vec, gi)
    push!(grad_ops, p)
end

# Line 220: Optimize angle for each generator
θi, costi = DBF.optimize_theta_expval(O, Gi, ψ, verbose=0)

# Line 235: Evolve operator
evolve!(O, Gi, θi)
```

**Optimization Function:** `optimize_theta_expval` (line 11)
- Computes: `E(θ) = cos²(θ/2)⟨ψ|H|ψ⟩ + sin²(θ/2)⟨Gψ|H|Gψ⟩ - 2i cos(θ/2)sin(θ/2)⟨ψ|H|Gψ⟩`
- Uses Optim.jl to find minimum in range [0, 2π]

---

### Canonical DBF Implementation

**Location:** `src/groundstate.jl:411` - function `groundstate_diffeq`

**Logic Flow:**
```julia
1. Create reference projector: S = create_0_projector(N, n_body)
2. Repeat until converged:
   a. Compute commutator: pool = S*O - O*S  # [S, H]
   b. Truncate by coefficient
   c. Extract gradient directly from commutator coefficients:
      For each term (p, c) in pool:
         θ_i = -Im(c)  # Direct gradient (no state evaluation)
   d. For all gradient terms:
      - Use: θ_i ← θ_i × stepsize
      - Evolve: H ← exp(iθ_i G_i/2) H exp(-iθ_i G_i/2)
      - Truncate H
   e. Check convergence: ||∇|| < threshold
```

**Key Code Sections:**
```julia
# Line 435: Create projector (same as Driven)
S = create_0_projector(N, n_body)

# Line 457: Compute commutator directly
pool = S*O - O*S

# Line 473-482: Extract angles from imaginary coefficients
for (p, c) in pool
    push!(grad_vec, -imag(c))  # θ = -Im(c)
    push!(grad_ops, p)
end

# Line 519: Evolve with scaled angle
O = evolve(O, gi, θi*stepsize)  # Note: θi*stepsize, not optimized θ
```

**Key Difference:**
- **No angle optimization** - uses commutator coefficients directly
- **Fixed stepsize** parameter controls rotation magnitude
- **Faster per iteration** but may need more iterations

---

### Diagonalization DBF Implementation

**Location:** `src/diagonalization.jl:12` - function `dbf_diag`

**Logic Flow:**
```julia
1. Repeat until converged:
   a. Extract diagonal: S = diag(O)  # Changes each iteration
   b. Compute top commutator terms: SO_OS = max_of_commutator2(S, O, n_top)
   c. Find largest: (coeff, G) = findmax(|·|, SO_OS)
   d. Optimize angle: θ_i = argmax_θ ||diag(H(θ))||²
   e. Evolve: H ← exp(iθ_i G/2) H exp(-iθ_i G/2)
   f. Truncate H
   g. Check: diagonal norm increasing, no repeated generators
```

**Key Code Sections:**
```julia
# Line 29: Extract diagonal (changes each iteration!)
S = diag(O)

# Line 45: Find top commutator terms (efficiency optimization)
SO_OS = max_of_commutator2(S, O, n_top=search_n_top)

# Line 52: Select largest commutator
coeff, G = findmax(v -> abs(v), SO_OS)

# Line 56: Optimize for diagonalization (different cost!)
θi, costi = DBF.optimize_theta_diagonalization(O, G, verbose=0)

# Line 62: Evolve
O = evolve(O, G, θi)
```

**Optimization Function:** `optimize_theta_diagonalization` (line 147)
- Computes: `C(θ) = ||diag(H(θ))||²`
- Analytical formula using inner products of diagonal parts
- Maximizes diagonal content (not minimizes energy)

**Cost Function Details:**
```julia
# Line 154-159: Compute inner products
A = inner_product(Od, Od)        # ⟨diag(H)|diag(H)⟩
B = inner_product(GOGd, GOGd)    # ⟨diag(GHG)|diag(GHG)⟩
C = inner_product(comd, comd)    # ⟨diag([G,H])|diag([G,H])⟩
D = inner_product(Od, GOGd)      # ⟨diag(H)|diag(GHG)⟩
E = inner_product(Od, comd)      # ⟨diag(H)|diag([G,H])⟩
F = inner_product(GOGd, comd)    # ⟨diag(GHG)|diag([G,H])⟩

# Line 208-210: Cost function
c = cos(θ/2)
s = sin(θ/2)
cost = c^4*A + s^4*B + s^2*c^2*(C+2D) + 2im*c^3*s*E + 2im*c*s^3*F
```

---

### Disentanglement DBF Implementation

**Location:** `src/disentangle.jl:98` - function `dbf_disentangle`

**Logic Flow:**
```julia
1. Define split index M (separates P and Q spaces)
2. Repeat until converged:
   a. Project to P-space: source = p_space(O, M)  # Changes each iteration
   b. Compute commutator: com = O*source - source*O
   c. Find largest: (coeff, G) = findmax(|·|, com)
   d. Optimize angle: θ_i = argmin_θ ||q_space(H(θ), M)||²
   e. Evolve: H ← exp(iθ_i G/2) H exp(-iθ_i G/2)
   f. Truncate H
   g. Check convergence: ||Q(H)|| < threshold
```

**Key Code Sections:**
```julia
# Line 121: Project to P-space (changes each iteration!)
source = p_space(O, M)

# Line 122: Compute commutator
com = O*source - source*O

# Line 128: Find largest commutator
coeff, G = findmax(v -> abs(v), com)

# Line 129: Optimize for disentanglement (different cost!)
θi, costi = optimize_theta_disentangle(O, G, M, verbose=0)

# Line 130: Evolve
O = evolve(O, G, θi)

# Line 133: Check Q-space norm
norm_new = norm(q_space(O, M))
```

**Space Projection Functions:**

`p_space` (line 2):
```julia
function p_space(O::PauliSum{N,T}, split_idx) where {N,T}
    # Returns operators acting only on qubits ≤ split_idx
    # Condition: (p.x|p.z) < 2^M OR acts entirely on lower qubits
end
```

`q_space` (line 13):
```julia
function q_space(O::PauliSum{N,T}, split_idx) where {N,T}
    # Returns operators coupling across the split
    # Condition: (p.x|p.z) >= 2^M AND acts across boundary
end
```

**Optimization Function:** `optimize_theta_disentangle` (line 24)
- Computes: `C(θ) = ||P(H(θ))||²` 
- Similar form to diagonalization but projects onto P-space instead of diagonal
- Maximizes P-space norm (equivalently minimizes Q-space coupling across the cut)

---

## Convergence Criteria Comparison

### Driven DBF
```julia
# Converges when gradient norm is small
if norm(grad_vec) < conv_thresh
    break
end
```
- **Measures:** `||∇|| = ||⟨[P,H]⟩||` (norm of gradient vector)
- **Typical threshold:** 1e-3

### Canonical DBF
```julia
# Converges when commutator norm is small
norm_new = norm(grad_vec)
if norm_new < conv_thresh
    break
end
```
- **Measures:** `||[P,H]||` (commutator norm)
- **Typical threshold:** 1e-3

### Diagonalization DBF
```julia
# Stops when diagonal norm stops increasing
if norm_new > norm_old
    @warn "Norm increased?"
    break
end
```
- **Measures:** Whether `||diag(H)||` is increasing
- Also stops if generator repeats (stuck in cycle)

### Disentanglement DBF
```julia
# Converges when Q-space norm is small
norm_new = norm(q_space(O, M))
if norm_new < conv_thresh
    break
end
```
- **Measures:** `||Q(H)||` (Q-space norm)
- **Typical threshold:** 1e-3

---

## When to Use Each Variant

### Use Driven DBF (`dbf_groundstate`) when:
- ✅ Finding ground state energy in a specific reference state
- ✅ You want high accuracy per rotation (willing to pay optimization cost)
- ✅ You have a good reference state (Hartree-Fock, product state, etc.)
- ✅ Energy minimization is the primary goal
- ✅ You want to track energy, variance, PT2 corrections

**Example Use Cases:**
- Finding ground state of molecular Hamiltonians
- Variational quantum eigensolver (VQE) type algorithms
- When reference state is chemically or physically motivated

### Use Canonical DBF (`groundstate_diffeq`) when:
- ✅ You want faster iterations (no per-generator optimization)
- ✅ Willing to take more total rotations
- ✅ Exploring stepsize/discretization effects
- ✅ Research into differential equation formulation

**Example Use Cases:**
- Rapid prototyping
- Large systems where optimization cost is prohibitive
- Studying continuous flow limit

### Use Diagonalization DBF (`dbf_diag`) when:
- ✅ Need eigenvalue spectrum (not just ground state)
- ✅ Want to see full energy landscape
- ✅ No specific reference state required
- ✅ Interested in all eigenstates

**Example Use Cases:**
- Finding all eigenvalues of small Hamiltonians
- Comparing with exact diagonalization
- Studying spectral properties
- Pre-processing before other algorithms

### Use Disentanglement DBF (`dbf_disentangle`) when:
- ✅ Want to separate quantum system into subsystems
- ✅ Reducing entanglement across a boundary
- ✅ Preparing for tensor network methods
- ✅ Studying subsystem structure

**Example Use Cases:**
- Preparing for DMRG (Density Matrix Renormalization Group)
- Studying entanglement structure
- Partitioning chemical systems into fragments
- Quantum circuit optimization (reducing two-qubit gates across regions)

---

## Computational Cost Comparison

| Variant | Cost per Iteration | Iterations Needed | Total Cost | Notes |
|---------|-------------------|-------------------|------------|-------|
| **Driven** | High (angle opt) | Medium | Medium-High | Most accurate per rotation |
| **Canonical** | Low (no opt) | High | Medium | Fastest per iteration |
| **Diagonalization** | High (angle opt) | Medium-High | High | Complete spectrum |
| **Disentanglement** | High (angle opt) | Medium | Medium-High | Depends on entanglement |

**Per-iteration breakdown:**

Driven DBF:
- Commutator: `O(|P| × |H|)`
- Gradient vector: `O(|G| × |basis|)` with matvec
- Per-generator optimization: `O(100 × |H|)` (100 function evaluations typical)
- Evolution: `O(|G| × |H|)`

Canonical DBF:
- Commutator: `O(|P| × |H|)`
- No optimization!
- Evolution: `O(|G| × |H|)` for all generators at once

Diagonalization DBF:
- Diagonal extraction: `O(|H|)`
- Top-k commutator: `O(k × k)` where `k = n_top`
- Per-generator optimization: `O(100 × |H|)`
- Evolution: `O(1 × |H|)` (only one generator per iteration)

Disentanglement DBF:
- P-space projection: `O(|H|)`
- Commutator: `O(|H_P| × |H|)`
- Per-generator optimization: `O(100 × |H|)`
- Evolution: `O(1 × |H|)`

---

## Output Format Comparison

### Driven DBF Returns:
```julia
Dict(
    "hamiltonian" => H_final,           # Evolved Hamiltonian
    "state" => ψ,                       # Reference state
    "energies" => Vector,               # Energy at each rotation
    "variances" => Vector,              # Variance at each rotation
    "accumulated_error" => Vector,       # Truncation error
    "generators" => Vector{PauliBasis}, # Generators used
    "angles" => Vector{Float64},        # Rotation angles
    "energies_per_grad" => Vector,      # Energy per gradient iteration
    "pt2_per_grad" => Vector,           # PT2 correction per iteration
    ...
)
```

### Canonical DBF Returns:
```julia
(O_final, generators, angles)
```
- Simpler output
- No detailed tracking

### Diagonalization DBF Returns:
```julia
(O_diag, generators, angles)
```
- `O_diag`: Diagonalized operator
- Eigenvalues can be extracted from diagonal terms

### Disentanglement DBF Returns:
```julia
(O_disentangled, generators, angles)
```
- `O_disentangled`: Block-diagonal operator
- P-space and Q-space are decoupled

---

## Example Comparison

Let's compare all four on the same Hamiltonian:

```julia
using DBF
using PauliOperators

N = 6
H = heisenberg_1D(N, 1.0, 1.0, 1.0)
ψ = Ket{N}(0)

println("Original H: $(length(H)) terms\n")

# 1. Driven DBF
result1 = dbf_groundstate(H, ψ, n_body=2, max_iter=20, verbose=0)
println("Driven DBF:")
println("  Final energy: $(result1["energies"][end])")
println("  Rotations: $(length(result1["angles"]))")
println("  Final terms: $(length(result1["hamiltonian"]))")

# 2. Canonical DBF  
H2, gen2, ang2 = groundstate_diffeq(H, ψ, n_body=2, max_iter=20, 
                                    stepsize=0.01, verbose=0)
E2 = expectation_value(H2, ψ)
println("\nCanonical DBF:")
println("  Final energy: $E2")
println("  Rotations: $(length(ang2))")
println("  Final terms: $(length(H2))")

# 3. Diagonalization DBF
H3, gen3, ang3 = dbf_diag(H, max_iter=20, verbose=0)
eigenvals = sort([c for (p,c) in diag(H3)], rev=true)
println("\nDiagonalization DBF:")
println("  Ground state energy: $(eigenvals[end])")
println("  Rotations: $(length(ang3))")
println("  Diagonal terms: $(length(diag(H3)))")
println("  Off-diagonal terms: $(length(offdiag(H3)))")

# 4. Disentanglement DBF
M = N ÷ 2  # Split in middle
H4, gen4, ang4 = dbf_disentangle(H, M, max_iter=20, verbose=0)
println("\nDisentanglement DBF:")
println("  Rotations: $(length(ang4))")
println("  P-space norm: $(norm(p_space(H4, M)))")
println("  Q-space norm: $(norm(q_space(H4, M)))")
```

**Expected Output Pattern:**
- **Driven:** Lowest energy, moderate rotations, tracks multiple metrics
- **Canonical:** Similar energy, more rotations, faster per iteration
- **Diagonalization:** Full spectrum, many rotations, most terms
- **Disentanglement:** System separated, Q-space small, different goal

---

## Summary

The four DBF variants share the same core evolution mechanism (`evolve` function) but differ fundamentally in:

1. **What they optimize:** Energy (Driven/Canonical), Diagonal norm (Diagonalization), Q-space norm (Disentanglement)

2. **How they compute gradients:** State-based (Driven), Direct commutator (Canonical), Diagonal-based (Diagonalization), P-space projection (Disentanglement)

3. **How they choose angles:** Optimized per generator (Driven/Diag/Disent), Fixed stepsize (Canonical)

4. **What drives evolution:** Fixed projector (Driven/Canonical), Changing diagonal (Diagonalization), Changing P-space (Disentanglement)

5. **What they output:** Energy trajectory (Driven), Simple final state (Canonical), Eigenspectrum (Diagonalization), Block structure (Disentanglement)

Choose the variant based on your computational goal!
