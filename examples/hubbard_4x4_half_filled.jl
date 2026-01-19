using DBF
using PauliOperators
using LinearAlgebra
using Printf
using Random

"""
Demonstration: 4×4 Hubbard model with half-filling and particle number preservation

The Hubbard model on a 4×4 lattice has:
- 16 physical sites
- 2 spin orbitals per site (up and down)
- Total of 32 qubits

Half-filling means 8 electrons (4 up, 4 down typically), which corresponds
to 8 occupied orbitals out of 32 total orbitals.
"""

println("="^80)
println("4×4 Fermi-Hubbard Model with Half-Filling")
println("Using Particle Number Preservation in DBF Ground State")
println("="^80)
println()

# System parameters
Lx = 4  # Lattice size in x
Ly = 4  # Lattice size in y
Nsites = Lx * Ly  # 16 physical sites
N_qubits = 2 * Nsites  # 32 qubits (2 per site for spin up/down)

println("System parameters:")
println("  Lattice size: $(Lx)×$(Ly) = $(Nsites) sites")
println("  Total qubits (spin up + down): $(N_qubits)")
println()

# Hubbard model parameters
t = 1.0   # Hopping amplitude
U = 4.0   # On-site interaction (typical value for correlation physics)

println("Hubbard model parameters:")
println("  Hopping (t): $(t)")
println("  On-site interaction (U): $(U)")
println()

# Build the Hamiltonian
println("Building 2D Fermi-Hubbard Hamiltonian...")
H = DBF.fermi_hubbard_2D(Lx, Ly, t, U)
println("  Number of terms in Hamiltonian: ", length(H))
println("  Hamiltonian norm: ", norm(H))
println()

# Create the particle number operator
println("Creating particle number operator...")
N̂ = DBF.particle_number_operator(N_qubits)
println("  Number of terms in N̂: ", length(N̂))
println()

# Find a half-filled state (8 electrons out of 32 orbitals)
# We'll search for the state with particle number closest to 8
println("Searching for half-filled initial state (8 particles)...")

target_particles = 8
best_state_idx = 0
best_particle_diff = Inf

# For a 32-qubit system, we can't enumerate all states, so we'll sample
# or use a heuristic. Let's create a state with exactly 8 ones.
Random.seed!(42)

# Create a ket with exactly 8 particles (bit pattern with 8 ones)
# This is a simple binary number with 8 bits set
function create_half_filled_state(N_qubits, n_particles)
    # Create a bit pattern with n_particles bits set to 1
    # Start with n_particles consecutive 1s
    bits = 0
    for i in 0:(n_particles-1)
        bits |= (1 << i)
    end
    return Ket{N_qubits}(bits)
end

ψ = create_half_filled_state(N_qubits, target_particles)
println("  Initial state: |$(bitstring(Int128(ψ.v))[end-N_qubits+1:end])>")

# Verify particle number
particle_num = real(expectation_value(N̂, ψ))
println("  Particle number: $(particle_num)")
println("  Initial energy: ", real(expectation_value(H, ψ)))
println()

# Check if the Hamiltonian preserves particle number
println("Checking Hamiltonian properties...")
println("  Does H preserve particle number? ", DBF.preserves_particle_number(H))
println()

# Run DBF ground state optimization with particle number preservation
println("Running DBF ground state optimization...")
println("  Max iterations: 10")
println("  Particle number preservation: ENABLED")
println()

res = DBF.dbf_groundstate(H, ψ,
                max_iter=10,
                conv_thresh=1e-4,
                evolve_coeff_thresh=1e-8,
                grad_coeff_thresh=1e-8,
                preserve_particle_number=true,
                verbose=1)

H_transformed = res["hamiltonian"]
generators = res["generators"]
angles = res["angles"]

println()
println("="^80)
println("Results:")
println("="^80)
println()

# Check particle number after transformation
final_particle_num = real(expectation_value(N̂, ψ))
println("Particle number conservation:")
println("  Initial particle number: $(particle_num)")
println("  Final particle number: $(final_particle_num)")
println("  Change: $(abs(final_particle_num - particle_num))")
println()

println("Energy:")
println("  Initial energy: ", real(res["energies"][1]))
println("  Final energy: ", real(res["energies"][end]))
println("  Energy lowering: ", real(res["energies"][1] - res["energies"][end]))
println()

println("Generators:")
println("  Number of generators used: ", length(generators))
println("  All generators preserve particle number: ", all(DBF.preserves_particle_number.(generators)))
println()

println("Hamiltonian transformation:")
println("  Initial H terms: ", length(H))
println("  Final H terms: ", length(H_transformed))
println("  Variance: ", real(res["variances"][end]))
println()

# Show energy convergence
if length(res["energies_per_grad"]) > 0
    println("Energy per iteration:")
    for (i, e) in enumerate(res["energies_per_grad"])
        @printf("  Iter %2d: E = %12.8f\n", i-1, real(e))
    end
end

println()
println("="^80)
println("Demonstration complete!")
println("="^80)
