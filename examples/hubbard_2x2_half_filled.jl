using DBF
using PauliOperators
using LinearAlgebra
using Printf
using Random

"""
Demonstration: 2×2 Hubbard model with half-filling and particle number preservation

For demonstration purposes, we use a smaller 2×2 system which is more tractable:
- 4 physical sites
- 2 spin orbitals per site (up and down)
- Total of 8 qubits

Half-filling means 2 electrons (1 up, 1 down or 2 up, etc.), which corresponds
to 2 occupied orbitals out of 8 total orbitals.
"""

println("="^80)
println("2×2 Fermi-Hubbard Model with Half-Filling")
println("Using Particle Number Preservation in DBF Ground State")
println("="^80)
println()

# System parameters
Lx = 2  # Lattice size in x
Ly = 2  # Lattice size in y
Nsites = Lx * Ly  # 4 physical sites
N_qubits = 2 * Nsites  # 8 qubits (2 per site for spin up/down)

println("System parameters:")
println("  Lattice size: $(Lx)×$(Ly) = $(Nsites) sites")
println("  Total qubits (spin up + down): $(N_qubits)")
println()

# Hubbard model parameters
t = 1.0   # Hopping amplitude
U = 4.0   # On-site interaction

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

# For 8 qubits, we can search through states to find a good initial state
# with the target particle number
target_particles = 2
println("Searching for initial state with $(target_particles) particles...")

# Search through all states with exactly target_particles bits set
function find_best_initial_state(H, N̂, N_qubits, target_particles)
    best_energy = Inf
    best_state = nothing
    
    for state_idx in 0:(2^N_qubits - 1)
        if count_ones(state_idx) == target_particles
            ψ_test = Ket{N_qubits}(state_idx)
            particle_num = real(expectation_value(N̂, ψ_test))
            
            # Check if this state has the right particle number
            if abs(particle_num - target_particles) < 0.1
                energy = real(expectation_value(H, ψ_test))
                if energy < best_energy
                    best_energy = energy
                    best_state = ψ_test
                end
            end
        end
    end
    return best_state
end

ψ = find_best_initial_state(H, N̂, N_qubits, target_particles)
println("  Best initial state found: |$(bitstring(Int128(ψ.v))[end-N_qubits+1:end])>")

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
println("  Max iterations: 20")
println("  Particle number preservation: ENABLED")
println()

res = DBF.dbf_groundstate(H, ψ,
                max_iter=20,
                conv_thresh=1e-5,
                evolve_coeff_thresh=1e-10,
                grad_coeff_thresh=1e-10,
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
if length(generators) > 0
    println("  All generators preserve particle number: ", all(DBF.preserves_particle_number.(generators)))
    println("  First few generators:")
    for (i, gen) in enumerate(generators[1:min(5, length(generators))])
        println("    $(i): $(string(gen))")
    end
end
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
