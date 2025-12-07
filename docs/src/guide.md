# User Guide and Examples

## Getting Started

### Installation

```julia
using Pkg
Pkg.add("DBF")
```

Or install from the repository:
```julia
Pkg.add(url="https://github.com/arnab82/DBF.jl")
```

### Basic Setup

```julia
using DBF
using PauliOperators
```

## Core Concepts

### Pauli Operators

Pauli operators are the building blocks of quantum mechanics in the qubit representation:

```julia
N = 4  # Number of qubits

# Identity
I = Pauli(N)

# Single-qubit Paulis
X1 = Pauli(N, X=[1])      # X on qubit 1
Z2 = Pauli(N, Z=[2])      # Z on qubit 2
Y3 = Pauli(N, Y=[3])      # Y on qubit 3

# Multi-qubit Paulis
X1X2 = Pauli(N, X=[1,2])  # X₁ ⊗ X₂
Z1Z3 = Pauli(N, Z=[1,3])  # Z₁ ⊗ Z₃
Y1X2Z3 = Pauli(N, Y=[1], X=[2], Z=[3])  # Y₁ ⊗ X₂ ⊗ Z₃
```

### Pauli Sums (Hamiltonians)

Operators are represented as weighted sums of Pauli strings:

```julia
H = PauliSum(N, Float64)

# Add terms
H += 0.5 * Pauli(N, X=[1,2])
H += -1.0 * Pauli(N, Z=[1])
H += 0.25 * Pauli(N, Y=[2], Z=[3])

# Or build from scratch
H = 0.5 * Pauli(N, X=[1,2]) + 
    -1.0 * Pauli(N, Z=[1]) + 
    0.25 * Pauli(N, Y=[2], Z=[3])
```

### Quantum States

States are represented as linear combinations of computational basis states:

```julia
# Computational basis states
ψ0 = Ket{N}(0)      # |0000⟩
ψ1 = Ket{N}(1)      # |0001⟩
ψ5 = Ket{N}(5)      # |0101⟩
ψall = Ket{N}(15)   # |1111⟩

# Superposition states
ψ = KetSum(N, ComplexF64)
ψ[Ket{N}(0)] = 1/√2
ψ[Ket{N}(15)] = 1/√2
# This is (|0000⟩ + |1111⟩)/√2
```

### Expectation Values

```julia
# Energy expectation value
E = expectation_value(H, ψ)

# Variance
V = variance(H, ψ)
```

## Example 1: Basic Pauli Evolution

Evolving a Hamiltonian under a unitary transformation:

```julia
using DBF
using PauliOperators

N = 4
H = PauliSum(N, Float64)
H += 1.0 * Pauli(N, X=[1,2])
H += 0.5 * Pauli(N, Z=[1])
H += -0.3 * Pauli(N, Y=[2,3])

# Generator for evolution
G = PauliBasis(Pauli(N, Z=[1,2]))

# Rotation angle
θ = π/4

# Evolve: H' = exp(iθG/2) H exp(-iθG/2)
H_evolved = evolve(H, G, θ)

println("Original H:")
display(H)
println("\nEvolved H:")
display(H_evolved)
```

## Example 2: Ground State Finding

Finding the ground state energy of a Heisenberg chain:

```julia
using DBF
using PauliOperators

# 1D Heisenberg chain with 6 qubits
N = 6
J = 1.0
H = heisenberg_1D(N, J, J, J)

println("Hamiltonian has $(length(H)) terms")

# Reference state (all spins down)
ψ = Ket{N}(0)

# Initial energy
E0 = expectation_value(H, ψ)
println("Initial energy: $E0")

# Run DBF ground state algorithm
result = dbf_groundstate(H, ψ,
    n_body = 2,                    # 2-body projector approximation
    max_iter = 100,                # Maximum iterations
    evolve_coeff_thresh = 1e-8,    # Truncation threshold
    grad_coeff_thresh = 1e-6,      # Gradient threshold
    evolve_weight_thresh = N,      # Max Pauli weight
    max_rots_per_grad = 50,        # Max rotations per gradient
    verbose = 1)                   # Print progress

# Extract results
E_final = result["energies"][end]
H_final = result["hamiltonian"]
n_rotations = length(result["angles"])

println("\nFinal energy: $E_final")
println("Number of rotations: $n_rotations")
println("Energy lowering: $(E0 - E_final)")

# Check PT2 correction
e0, e2 = pt2(H_final, ψ)
println("PT2 correction: $e2")
println("Total energy estimate: $(e0 + e2)")
```

## Example 3: Hamiltonian Diagonalization

Diagonalizing a spin Hamiltonian:

```julia
using DBF
using PauliOperators

N = 6
H = heisenberg_1D(N, 1.0, 1.0, 1.0)

println("Initial operator:")
println("  Total terms: $(length(H))")
println("  Diagonal terms: $(length(diag(H)))")
println("  Off-diagonal terms: $(length(offdiag(H)))")
println("  Diagonal norm: $(norm(diag(H)))")
println("  Off-diagonal norm: $(norm(offdiag(H)))")

# Diagonalize
H_diag, generators, angles = dbf_diag(H,
    max_iter = 50,
    evolve_coeff_thresh = 1e-10,
    search_n_top = 100,
    verbose = 1)

println("\nFinal operator:")
println("  Total terms: $(length(H_diag))")
println("  Diagonal terms: $(length(diag(H_diag)))")
println("  Off-diagonal terms: $(length(offdiag(H_diag)))")
println("  Diagonal norm: $(norm(diag(H_diag)))")
println("  Off-diagonal norm: $(norm(offdiag(H_diag)))")

# The diagonal values are the eigenvalues
println("\nEigenvalue spectrum (sorted):")
eigenvalues = sort([c for (p,c) in diag(H_diag)], rev=true)
for (i, e) in enumerate(eigenvalues[1:min(10, end)])
    println("  λ[$i] = $e")
end
```

## Example 4: 2D Heisenberg Model

Working with 2D systems:

```julia
using DBF
using PauliOperators

# 3x3 lattice (9 qubits)
Nx, Ny = 3, 3
Jx, Jy, Jz = 1.0, 1.0, 1.0

# Create 2D Heisenberg Hamiltonian
H = heisenberg_2D(Nx, Ny, Jx, Jy, Jz, periodic=true)

println("2D Heisenberg ($Nx × $Ny):")
println("  Number of qubits: $(Nx*Ny)")
println("  Number of terms: $(length(H))")

# Reference state (Néel state)
ψ = Ket{Nx*Ny}(0b101010101)  # Alternating pattern

# Ground state search
result = dbf_groundstate(H, ψ,
    n_body = 2,
    max_iter = 50,
    evolve_coeff_thresh = 1e-6,
    verbose = 1)

E_final = result["energies"][end]
println("\nGround state energy: $E_final")
println("Energy per site: $(E_final / (Nx*Ny))")
```

## Example 5: Fermionic Systems (Hubbard Model)

Using the Fermi-Hubbard model:

```julia
using DBF
using PauliOperators

# 1D Hubbard chain with 4 sites (8 qubits for spin up/down)
L = 4
t = 1.0  # Hopping
U = 2.0  # On-site interaction

H = hubbard_model_1D(L, t, U)

println("Fermi-Hubbard 1D:")
println("  Sites: $L")
println("  Qubits: $(2*L)")
println("  Terms: $(length(H))")

# Half-filling: alternating up/down electrons
# |↑↓↑↓⟩ = |10101010⟩
ψ = Ket{2*L}(0b10101010)

E0 = expectation_value(H, ψ)
println("Initial energy: $E0")

# Run optimization
result = dbf_groundstate(H, ψ,
    n_body = 2,
    max_iter = 100,
    evolve_coeff_thresh = 1e-8,
    grad_coeff_thresh = 1e-6,
    verbose = 1)

println("\nFinal energy: $(result["energies"][end])")
```

## Example 6: ADAPT-VQE

Using ADAPT to build a variational ansatz:

```julia
using DBF
using PauliOperators

N = 6
H = heisenberg_1D(N, 1.0, 1.0, 1.0)

# Create operator pool (qubit excitations)
pool = Vector{PauliBasis{N}}()

# Add 1-body terms
for i in 1:N
    push!(pool, PauliBasis(Pauli(N, X=[i])))
    push!(pool, PauliBasis(Pauli(N, Y=[i])))
end

# Add 2-body terms
for i in 1:N
    for j in i+1:N
        push!(pool, PauliBasis(Pauli(N, Y=[i], X=[j])))
    end
end

println("Pool size: $(length(pool))")

# Reference state
ψ = Ket{N}(0)

# Run ADAPT
H_final, generators, angles = adapt(H, pool, ψ,
    max_iter = 20,
    grad_coeff_thresh = 1e-6,
    evolve_coeff_thresh = 1e-10,
    verbose = 1)

E_final = expectation_value(H_final, ψ)
println("\nFinal energy: $E_final")
println("Circuit depth: $(length(generators))")

println("\nCircuit structure:")
for (i, (g, θ)) in enumerate(zip(generators, angles))
    println("  $i: exp(-i θ/2 $(string(g))), θ = $(round(θ, digits=4))")
end
```

## Example 7: State Preparation

Preparing specific quantum states:

```julia
using DBF
using PauliOperators

N = 6

# Prepare 1D cluster state
generators, angles = get_1d_cluster_state_sequence(N)

println("1D Cluster State Preparation:")
println("  Number of gates: $(length(generators))")

# Start from |000000⟩
ψ = Ket{N}(0)
ψ_sum = KetSum(N, ComplexF64)
ψ_sum[ψ] = 1.0

# Apply gates in sequence
for (g, θ) in zip(generators, angles)
    ψ_sum = evolve(ψ_sum, g, θ)
end

println("Final state has $(length(ψ_sum)) basis states")

# Normalize
ψ_norm = norm(ψ_sum)
for (k, c) in ψ_sum
    ψ_sum[k] = c / ψ_norm
end

# Display state
println("\nCluster state composition:")
sorted_states = sort(collect(ψ_sum), by=x->abs(x[2]), rev=true)
for (k, c) in sorted_states[1:min(10, end)]
    bitstring = string(k.v, base=2, pad=N)
    println("  |$bitstring⟩: $(round(abs(c), digits=4))")
end
```

## Example 8: Truncation and Error Tracking

Understanding truncation effects:

```julia
using DBF
using PauliOperators

N = 8
H = heisenberg_1D(N, 1.0, 1.0, 1.0)
ψ = Ket{N}(0)

# Run with different truncation thresholds
thresholds = [1e-6, 1e-8, 1e-10, 1e-12]

println("Truncation threshold comparison:\n")
println("Threshold    Final E    Final Terms    Total Error")
println("="^60)

for thresh in thresholds
    result = dbf_groundstate(H, ψ,
        n_body = 2,
        max_iter = 20,
        evolve_coeff_thresh = thresh,
        grad_coeff_thresh = thresh * 100,
        verbose = 0)
    
    E_final = result["energies"][end]
    n_terms = length(result["hamiltonian"])
    total_error = result["accumulated_error"][end]
    
    @printf("%.0e    %10.6f    %8d    %12.8e\n", 
            thresh, E_final, n_terms, total_error)
end
```

## Example 9: Quantum Gates

Using quantum gates to manipulate operators and states:

```julia
using DBF
using PauliOperators

N = 4

# Start with simple operator
O = PauliSum(N, Float64)
O += 1.0 * Pauli(N, X=[1])

println("Initial operator: ", string(collect(keys(O))[1]))

# Apply Hadamard to qubit 1
O = hadamard(O, 1)
println("After H(1): ")
display(O)

# Apply CNOT from qubit 1 to qubit 2
O = cnot(O, 1, 2)
println("\nAfter CNOT(1,2): ")
display(O)

# Apply T gate to qubit 1
O = T_gate(O, 1)
println("\nAfter T(1): ")
display(O)

# For states:
ψ = KetSum(N, ComplexF64)
ψ[Ket{N}(0)] = 1.0

println("\n\nInitial state: |0000⟩")

# Create Bell pair
ψ = hadamard(ψ, 1)
ψ = cnot(ψ, 1, 2)

println("After H(1) and CNOT(1,2):")
for (k, c) in ψ
    if abs(c) > 1e-10
        bitstring = string(k.v, base=2, pad=N)
        println("  |$bitstring⟩: $(c)")
    end
end
# Should give (|0000⟩ + |1100⟩)/√2
```

## Example 10: Advanced - Subspace Methods

Using subspace methods for larger systems:

```julia
using DBF
using PauliOperators

N = 8
H = heisenberg_1D(N, 1.0, 1.0, 1.0)
ψ_ref = Ket{N}(0)

# FOIS-CI: First-Order Interacting Space CI
println("Running FOIS-CI...")

e0, e_vals, v_vecs, basis = fois_ci(H, ψ_ref,
    thresh = 1e-6,      # Basis truncation
    krylov_order = 2,   # H^2 |ψ⟩
    max_iter = 100,     # Eigensolver iterations
    tol = 1e-8,         # Convergence tolerance
    verbose = 1)

println("\nResults:")
println("  Reference energy: $e0")
println("  Ground state energy: $(e_vals[1])")
println("  Basis size: $(length(basis))")
println("  Correlation energy: $(e_vals[1] - e0)")

# CEPA: Coupled Electron Pair Approximation
println("\n\nRunning CEPA...")

e0, e_cepa, x, basis = cepa(H, ψ_ref,
    thresh = 1e-6,
    tol = 1e-8,
    verbose = 1)

println("\nResults:")
println("  Reference energy: $e0")
println("  CEPA energy: $(real(e_cepa))")
println("  Correlation energy: $(real(e_cepa - e0))")
```

## Tips and Best Practices

### Choosing Parameters

1. **n_body for Ground State:**
   - `n_body=1`: Fast, less accurate
   - `n_body=2`: Good balance (recommended)
   - `n_body≥3`: More accurate but expensive

2. **Truncation Thresholds:**
   - Start with `1e-6` for exploration
   - Use `1e-8` to `1e-12` for production
   - Monitor accumulated error

3. **Weight Limits:**
   - Start with `max_weight = N`
   - Reduce if operator grows too large
   - Use Majorana weight for fermions

4. **Convergence:**
   - `conv_thresh = 1e-3` for gradient norm
   - Increase `max_iter` if not converging
   - Check variance at convergence

### Performance Optimization

1. **Use Threading:**
   ```bash
   export JULIA_NUM_THREADS=4
   julia script.jl
   ```

2. **Checkpointing:**
   ```julia
   result = dbf_groundstate(H, ψ,
       checkfile = "checkpoint",
       ...)
   ```

3. **Monitor Progress:**
   ```julia
   verbose = 2  # Detailed output
   ```

### Debugging

1. **Check Operator Size:**
   ```julia
   println("Operator size: $(length(H))")
   println("Weight distribution: $(get_weight_counts(H))")
   ```

2. **Verify Energy Conservation:**
   ```julia
   E_before = expectation_value(H, ψ)
   H_evolved = evolve(H, G, θ)
   E_after = expectation_value(H_evolved, ψ)
   @assert abs(E_before - E_after) < 1e-10
   ```

3. **Check Hermiticity:**
   ```julia
   @assert inner_product(H, H) ≈ norm(H)^2
   ```

## Common Issues and Solutions

### Issue: Operator grows too large

**Solution:** Reduce truncation threshold or weight limit
```julia
evolve_coeff_thresh = 1e-6  # More aggressive truncation
evolve_weight_thresh = N-2   # Lower weight limit
```

### Issue: Not converging

**Solution:** 
- Increase `max_iter`
- Adjust `grad_coeff_thresh`
- Try different `n_body`
- Check initial state quality

### Issue: High accumulated error

**Solution:**
- Use tighter `evolve_coeff_thresh`
- Increase `max_weight`
- Use variance error tracking

## Further Examples

See the `test/` directory for more examples:
- `test/test_groundstate_dbf.jl` - Ground state examples
- `test/test_diag_dbf.jl` - Diagonalization examples
- `test/test_adapt.jl` - ADAPT examples
- `test/test_evolve.jl` - Evolution and gate examples

## Next Steps

1. Read the [Theory](theory.md) document for mathematical details
2. Explore the [Structure](structure.md) document for API reference
3. Try the examples with different Hamiltonians
4. Experiment with different parameters
5. Contribute your own examples!
