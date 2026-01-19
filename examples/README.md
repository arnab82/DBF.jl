# DBF.jl Examples

This directory contains examples demonstrating the use of DBF.jl for quantum many-body systems, with a focus on particle number preservation.

## Particle Number Preservation

The particle number preservation feature allows `dbf_groundstate` to restrict the search space to operators that preserve the total particle number in fermionic systems. This is particularly useful for:

- Ensuring physical constraints are maintained during optimization
- Reducing the effective search space for better efficiency
- Obtaining ground states within specific particle number sectors

## Examples

### 1. 1D Hubbard Model Comparison (`hubbard_1d_comparison.jl`)

Demonstrates the difference between running `dbf_groundstate` with and without particle number preservation on a 4-site 1D Hubbard model.

**Key features:**
- Compares optimization with and without particle number preservation
- Shows how generators are filtered when preservation is enabled
- Demonstrates that all generators preserve particle number when the flag is set

**Run:**
```bash
julia --project=. examples/hubbard_1d_comparison.jl
```

### 2. 2×2 Hubbard Model (`hubbard_2x2_half_filled.jl`)

Runs `dbf_groundstate` on a 2×2 square lattice Hubbard model with half-filling (2 electrons).

**System parameters:**
- 4 physical sites (2×2 lattice)
- 8 qubits total (2 spin orbitals per site)
- 2 electrons (half-filling)

**Run:**
```bash
julia --project=. examples/hubbard_2x2_half_filled.jl
```

### 3. 4×4 Hubbard Model (`hubbard_4x4_half_filled.jl`)

Demonstrates `dbf_groundstate` on a larger 4×4 square lattice Hubbard model with half-filling (8 electrons).

**System parameters:**
- 16 physical sites (4×4 lattice)
- 32 qubits total (2 spin orbitals per site)
- 8 electrons (half-filling)

**Run:**
```bash
julia --project=. examples/hubbard_4x4_half_filled.jl
```

## Understanding Particle Number Preservation

### What is the particle number operator?

For a system of N qubits, the particle number operator is:

```
N̂ = Σᵢ (I - Zᵢ)/2
```

This operator counts how many qubits are in the |1⟩ state (occupied orbitals).

### Which operators preserve particle number?

A Pauli operator preserves particle number if it commutes with N̂. For Pauli strings:

- **I and Z operators**: Always preserve particle number (they don't flip qubits)
- **X and Y operators**: Individual X or Y operators do NOT preserve particle number
- **Pairs of X/Y operators**: An even number of X/Y operators CAN preserve particle number

### Using particle number preservation in your code

```julia
using DBF
using PauliOperators

# Create a Hamiltonian
H = DBF.fermi_hubbard_2D(2, 2, 1.0, 4.0)

# Create initial state with specific particle number
ψ = Ket{8}(0b00000011)  # 2 particles

# Run DBF with particle number preservation
res = DBF.dbf_groundstate(H, ψ,
            max_iter=10,
            preserve_particle_number=true,  # Enable preservation
            verbose=1)

# All generators will preserve particle number
all_preserve = all(DBF.preserves_particle_number.(res["generators"]))
```

## Output Interpretation

When running these examples, you'll see:

1. **System parameters**: Lattice size, number of qubits, particle number
2. **Hamiltonian properties**: Number of terms, norm, particle number preservation
3. **Optimization progress**: Energy per iteration, variance, number of generators
4. **Results**:
   - Energy lowering achieved
   - Number of generators used
   - Verification that all generators preserve particle number
   - Final Hamiltonian statistics

## Additional Notes

- The Hubbard model Hamiltonians naturally preserve particle number (they commute with N̂)
- Initial states should have the desired particle number for physically meaningful results
- Particle number preservation is enforced by filtering gradient operators during optimization
- The feature adds minimal computational overhead (just filtering)
